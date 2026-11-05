
# PLATELAB - qPCR Data Analysis Website

import streamlit as st
import pandas as pd
from io import BytesIO
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

st.set_page_config(page_title="PlateLab", layout="wide")

st.title("PlateLab")
st.write("qPCR Experiment Design & Analysis")

# Upload data
uploaded_file = st.file_uploader(
    "Upload your qPCR results",
    type=["csv"]
)

if uploaded_file is not None:

    # STEP 1: Read and clean data
    raw_df = pd.read_csv(uploaded_file)

    raw_df.columns = raw_df.columns.str.strip().str.lower()

    raw_df = raw_df.rename(columns={
        "sample name": "Sample",
        "target name": "Gene",
        "ct": "Ct"
    })

    required = ["Sample", "Gene", "Ct"]

    if not all(c in raw_df.columns for c in required):
        st.error("Your CSV needs Sample, Gene, and Ct columns.")
        st.stop()

    raw_df["Sample"] = raw_df["Sample"].astype(str).str.strip()
    raw_df["Gene"] = raw_df["Gene"].astype(str).str.strip()
    raw_df["Ct"] = pd.to_numeric(raw_df["Ct"], errors="coerce")

    # STEP 2: Average technical replicates
    average_df = (
        raw_df.groupby(["Sample", "Gene"])
        .agg(
            Average_Ct=("Ct", "mean"),
            Valid_Replicates=("Ct", "count"),
            Total_Replicates=("Ct", "size"),
            Ct_SD=("Ct", "std")
        )
        .reset_index()
    )

    average_df = average_df.rename(columns={
        "Average_Ct": "Average Ct",
        "Valid_Replicates": "Valid Replicates",
        "Total_Replicates": "Total Replicates",
        "Ct_SD": "Ct SD"
    })

    # STEP 3: Analysis settings
    st.subheader("Analysis Settings")

    gene_options = sorted(raw_df["Gene"].dropna().unique())

    reference_gene = st.selectbox(
        "Select your reference gene",
        gene_options,
        index=gene_options.index("GAPDH")
        if "GAPDH" in gene_options else 0
    )

    sample_options = sorted(
        s for s in raw_df["Sample"].unique()
        if s.upper() != "NTC"
    )

    default_wt = [
        s for s in sample_options
        if s in ["WT 1", "WT 2", "WT 3"]
    ]

    reference_samples = st.multiselect(
        "Select samples for normalization",
        sample_options,
        default=default_wt
    )

    if not reference_samples:
        st.warning("Select at least one normalization sample.")
        st.stop()

    # STEP 4: Prepare analysis data
    analysis_df = average_df[
        average_df["Sample"].str.upper() != "NTC"
    ].copy()

    ref_df = analysis_df[
        analysis_df["Gene"] == reference_gene
    ][["Sample", "Gene", "Average Ct"]].copy()

    ref_df = ref_df.rename(columns={
        "Gene": "Reference Gene",
        "Average Ct": "Reference Ct"
    })

    target_df = analysis_df[
        analysis_df["Gene"] != reference_gene
    ].copy()

    target_df = target_df.rename(columns={
        "Gene": "Target Gene",
        "Average Ct": "Target Ct"
    })

    results = pd.merge(
        target_df,
        ref_df,
        on="Sample",
        how="left",
        validate="many_to_one"
    )

    # STEP 5: Delta Ct and relative expression
    results["Delta Ct"] = (
        results["Target Ct"] - results["Reference Ct"]
    )

    results["Relative Expression"] = (
        2 ** (-results["Delta Ct"])
    )

    # STEP 6: Normalize to selected samples
    control_df = results[
        results["Sample"].isin(reference_samples)
    ]

    control_average = (
        control_df.groupby("Target Gene")["Relative Expression"]
        .mean()
        .reset_index(name="Control Average")
    )

    normalized_df = pd.merge(
        results,
        control_average,
        on="Target Gene",
        how="left"
    )

    normalized_df["Normalized Expression"] = (
        normalized_df["Relative Expression"] /
        normalized_df["Control Average"]
    )

    # STEP 7: Clean summary
    summary_df = normalized_df[
        [
            "Sample",
            "Target Gene",
            "Delta Ct",
            "Relative Expression",
            "Normalized Expression"
        ]
    ].rename(columns={
        "Target Gene": "Gene"
    })

    # STEP 8: Display tabs
    tab1, tab2, tab3 = st.tabs([
        "Raw Data",
        "Average Ct & QC",
        "Results"
    ])

    with tab1:
        st.subheader("Raw qPCR Measurements")
        st.dataframe(raw_df, hide_index=True)

    with tab2:
        st.subheader("Average Ct & Quality Control")
        st.dataframe(
            average_df.round(4),
            hide_index=True
        )

        missing = raw_df["Ct"].isna().sum()
        st.write(f"Undetermined or missing Ct values: {missing}")

        st.info(
            "Review replicate variability and undetermined "
            "measurements before interpreting results."
        )

    with tab3:
        st.subheader("Normalized Gene Expression")

        st.write(f"Reference gene: {reference_gene}")
        st.write(
            "Normalization samples: " +
            ", ".join(reference_samples)
        )

        st.caption(
            "Method: 2^(-Delta Ct) divided by the "
            "arithmetic mean of 2^(-Delta Ct) "
            "for the selected reference samples."
        )

        st.dataframe(
            summary_df.round(6),
            hide_index=True
        )

        with st.expander("View complete calculations"):
            st.dataframe(
                normalized_df.round(6),
                hide_index=True
            )

        # STEP 9: Export Excel workbook
        buffer = BytesIO()

        with pd.ExcelWriter(
            buffer, engine="openpyxl"
        ) as writer:

           raw_df.to_excel(
    writer, sheet_name="Raw Ct Data", index=False
)

      average_df.to_excel(writer, sheet_name="Average Ct", index=False)

           normalized_df.to_excel(writer, sheet_name="Full Analysis", index=False)

            summary_df.to_excel(writer, sheet_name="Results Summary", index=False)

            # Style all sheets
            for sheet in writer.sheets.values():

                sheet.freeze_panes = "A2"
                sheet.auto_filter.ref = sheet.dimensions

                for cell in sheet[1]:
                    cell.font = Font(
                        bold=True, color="FFFFFF"
                    )
                    cell.fill = PatternFill(
                        fill_type="solid",
                        fgColor="24476B"
                    )

                for row in sheet.iter_rows():
                    for cell in row:
                        cell.alignment = Alignment(
                            horizontal="center",
                            vertical="center"
                        )

                        if isinstance(cell.value, float):
                            cell.number_format = "0.000000"

                for column in sheet.columns:
                    letter = get_column_letter(
                        column[0].column
                    )

                    max_length = max(
                        len(str(c.value))
                        if c.value is not None else 0
                        for c in column
                    )

                    sheet.column_dimensions[letter].width = min(
                        max(max_length + 3, 16), 32
                    )

        st.download_button(
            "Download Excel Results",
            data=buffer.getvalue(),
            file_name="PlateLab_Results.xlsx",
            mime=(
                "application/vnd.openxmlformats-officedocument."
                "spreadsheetml.sheet"
            )
        )
