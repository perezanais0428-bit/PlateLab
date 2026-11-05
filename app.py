# PLATELAB — qPCR Experiment Analysis Website
import streamlit as st
import pandas as pd
from io import BytesIO
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

st.set_page_config(page_title="PlateLab", page_icon="🧫", layout="wide")

# Visual styling only; does not affect the qPCR analysis.
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Archivo:wght@400;500;600;700&display=swap');
:root { --pl-pink: #eeb7d0; --pl-pink-hover: #e29cbc; }
html, body, .stApp, .stApp button, .stApp input,
.stApp textarea, .stApp label, .stApp p,
.stApp h1, .stApp h2, .stApp h3,
.stApp [data-testid="stMarkdownContainer"] {
  font-family: 'Archivo', sans-serif !important;
}
.stApp h1, .stApp h2, .stApp h3 { font-weight: 700 !important; }
/* Selected tabs and their underline */
.stTabs [data-baseweb="tab"] { font-family: 'Archivo', sans-serif !important; font-weight: 600; }
.stTabs [data-baseweb="tab"][aria-selected="true"] { color: var(--pl-pink) !important; }
.stTabs [data-baseweb="tab-highlight"] { background-color: var(--pl-pink) !important; }
.stTabs [data-baseweb="tab-border"] { opacity: .25; }
/* Selected normalization sample chips */
.stApp [data-baseweb="tag"] { background: var(--pl-pink) !important; border-radius: 9px !important; }
.stApp [data-baseweb="tag"] span { color: #392331 !important; font-family: 'Archivo', sans-serif !important; }
.stApp [data-baseweb="tag"] svg { color: #392331 !important; }
/* Input focus and dropdown outlines */
.stApp [data-baseweb="select"] > div:focus-within { border-color: var(--pl-pink) !important; box-shadow: 0 0 0 1px var(--pl-pink) !important; }
/* Download buttons */
.stApp [data-testid="stDownloadButton"] button { background-color: var(--pl-pink) !important; color: #392331 !important; border-color: var(--pl-pink) !important; border-radius: 10px !important; font-weight: 600 !important; }
.stApp [data-testid="stDownloadButton"] button:hover { background-color: var(--pl-pink-hover) !important; border-color: var(--pl-pink-hover) !important; }
.stApp [data-testid="stAlert"] { border-radius: 12px; }
</style>
""", unsafe_allow_html=True)

st.title("🧫 PlateLab")
st.write("qPCR Experiment Design & Analysis")
st.write("Upload your qPCR results to automatically calculate average Ct, Delta Ct, relative expression, and normalized gene expression.")

uploaded_file = st.file_uploader("Upload your qPCR CSV file", type=["csv"])

if uploaded_file is not None:
    raw_df = pd.read_csv(uploaded_file)
    raw_df.columns = raw_df.columns.str.strip().str.lower()
    raw_df = raw_df.rename(columns={
        "sample name": "Sample", "target name": "Gene", "ct": "Ct",
        "sample": "Sample", "gene": "Gene"
    })
    required = ["Sample", "Gene", "Ct"]
    missing_columns = [col for col in required if col not in raw_df.columns]
    if missing_columns:
        st.error("Missing columns: " + ", ".join(missing_columns))
        st.stop()

    raw_df = raw_df.dropna(subset=["Sample", "Gene"]).copy()
    raw_df["Sample"] = raw_df["Sample"].astype(str).str.strip()
    raw_df["Gene"] = raw_df["Gene"].astype(str).str.strip()
    raw_df["Ct"] = pd.to_numeric(raw_df["Ct"], errors="coerce")
    if raw_df.empty:
        st.error("No valid sample rows were found.")
        st.stop()

    average_df = (
        raw_df.groupby(["Sample", "Gene"], as_index=False)
        .agg(
            Average_Ct=("Ct", "mean"),
            Valid_Replicates=("Ct", "count"),
            Total_Replicates=("Ct", "size"),
            Ct_SD=("Ct", "std")
        )
    )
    average_df = average_df.rename(columns={
        "Average_Ct": "Average Ct", "Valid_Replicates": "Valid Replicates",
        "Total_Replicates": "Total Replicates", "Ct_SD": "Ct SD"
    })

    st.header("Analysis Settings")
    genes = sorted(raw_df["Gene"].unique())
    default_gene = genes.index("GAPDH") if "GAPDH" in genes else 0
    reference_gene = st.selectbox("Select Reference Gene", options=genes, index=default_gene)
    samples = sorted(sample for sample in raw_df["Sample"].unique() if sample.upper() != "NTC")
    wt_defaults = [sample for sample in samples if sample in ["WT 1", "WT 2", "WT 3"]]
    normalization_samples = st.multiselect(
        "Select Samples for Normalization", options=samples, default=wt_defaults
    )
    if not normalization_samples:
        st.warning("Please select at least one normalization sample.")
        st.stop()

    analysis_df = average_df[average_df["Sample"].str.upper() != "NTC"].copy()
    reference_df = analysis_df[analysis_df["Gene"] == reference_gene][
        ["Sample", "Gene", "Average Ct"]
    ].copy()
    reference_df = reference_df.rename(columns={
        "Gene": "Reference Gene", "Average Ct": "Reference Ct"
    })
    target_df = analysis_df[analysis_df["Gene"] != reference_gene].copy()
    target_df = target_df.rename(columns={
        "Gene": "Target Gene", "Average Ct": "Target Ct"
    })
    if target_df.empty:
        st.error("No target genes were found. Your file needs at least one gene besides the reference gene.")
        st.stop()
    results = pd.merge(target_df, reference_df, on="Sample", how="left", validate="many_to_one")
    results["Delta Ct"] = results["Target Ct"] - results["Reference Ct"]
    results["Relative Expression"] = 2 ** (-results["Delta Ct"])

    control_df = results[results["Sample"].isin(normalization_samples)]
    control_average = control_df.groupby("Target Gene", as_index=False)["Relative Expression"].mean()
    control_average = control_average.rename(columns={"Relative Expression": "Control Average"})
    normalized_df = pd.merge(results, control_average, on="Target Gene", how="left", validate="many_to_one")
    normalized_df["Normalized Expression"] = (
        normalized_df["Relative Expression"] / normalized_df["Control Average"]
    )
    invalid_control = (
        normalized_df["Control Average"].isna() |
        (normalized_df["Control Average"] <= 0)
    )
    normalized_df.loc[invalid_control, "Normalized Expression"] = float("nan")

    summary_df = normalized_df[[
        "Sample", "Target Gene", "Delta Ct", "Relative Expression", "Normalized Expression"
    ]].copy()
    summary_df = summary_df.rename(columns={
        "Target Gene": "Gene", "Normalized Expression": "Normalized Fold"
    })

    tab1, tab2, tab3 = st.tabs(["Raw Data", "Average Ct & QC", "Results"])
    with tab1:
        st.subheader("Raw qPCR Data")
        st.write(f"Total wells imported: {len(raw_df)}")
        st.dataframe(raw_df, hide_index=True, use_container_width=True)

    with tab2:
        st.subheader("Average Ct Results")
        st.dataframe(average_df.round(4), hide_index=True, use_container_width=True)
        missing_ct = raw_df["Ct"].isna().sum()
        st.write(f"Undetermined or missing Ct values: {missing_ct}")
        # QC review threshold is configurable in a future version.
        qc_flags = average_df[
            (average_df["Ct SD"] > 0.5) |
            (average_df["Valid Replicates"] < average_df["Total Replicates"])
        ]
        if not qc_flags.empty:
            st.warning("Some measurements need QC review.")
            st.dataframe(qc_flags.round(4), hide_index=True, use_container_width=True)
        else:
            st.success("No replicate variability or missing-value flags detected.")
        st.caption("QC flags are for review only. No Ct values are automatically removed. A passed check does not validate an assay.")

    with tab3:
        st.subheader("Gene Expression Results")
        st.write(f"Reference Gene: {reference_gene}")
        st.write("Normalization Samples: " + ", ".join(normalization_samples))
        st.caption(
            "Normalization method: each sample's 2^(-Delta Ct) divided by the arithmetic mean "
            "2^(-Delta Ct) of the selected reference samples, separately for each gene."
        )
        st.dataframe(summary_df.round(6), hide_index=True, use_container_width=True)
        with st.expander("View Complete Calculations"):
            st.dataframe(normalized_df.round(6), hide_index=True, use_container_width=True)
        if invalid_control.any():
            st.warning("Some genes have no valid normalization average. Review your selected reference samples and missing Ct values.")

        buffer = BytesIO()
        with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
            raw_df.to_excel(writer, sheet_name="Raw Ct Data", index=False)
            average_df.to_excel(writer, sheet_name="Average Ct", index=False)
            normalized_df.to_excel(writer, sheet_name="Full Analysis", index=False)
            summary_df.to_excel(writer, sheet_name="Results Summary", index=False)
            for sheet in writer.sheets.values():
                sheet.freeze_panes = "A2"
                sheet.auto_filter.ref = sheet.dimensions
                for cell in sheet[1]:
                    cell.font = Font(bold=True, color="FFFFFF")
                    cell.fill = PatternFill(fill_type="solid", fgColor="24476B")
                for row in sheet.iter_rows():
                    for cell in row:
                        cell.alignment = Alignment(horizontal="center", vertical="center")
                        if isinstance(cell.value, float):
                            cell.number_format = "0.000000"
                for column in sheet.columns:
                    letter = get_column_letter(column[0].column)
                    max_length = max(len(str(cell.value)) if cell.value is not None else 0 for cell in column)
                    sheet.column_dimensions[letter].width = min(max(max_length + 3, 16), 32)

        st.download_button(
            label="Download Excel Results",
            data=buffer.getvalue(),
            file_name="PlateLab_Results.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
