
# ==========================================
# PLATELAB
# qPCR Experiment Analysis Website
# ==========================================

import streamlit as st
import pandas as pd
from io import BytesIO
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


# ==========================================
# 1. WEBSITE SETUP
# ==========================================

st.set_page_config(
    page_title="PlateLab",
    page_icon="🧫",
    layout="wide"
)

st.title("🧫 PlateLab")
st.write("qPCR Experiment Design & Analysis")

st.write(
    "Upload your qPCR results to automatically "
    "calculate average Ct, Delta Ct, relative "
    "expression, and normalized gene expression."
)


# ==========================================
# 2. UPLOAD DATA
# ==========================================

uploaded_file = st.file_uploader(
    "Upload your qPCR CSV file",
    type=["csv"]
)

if uploaded_file is not None:

    # Read the uploaded CSV
    raw_df = pd.read_csv(uploaded_file)

    # Standardize column names
    raw_df.columns = (
        raw_df.columns.str.strip().str.lower()
    )

    raw_df = raw_df.rename(columns={
        "sample name": "Sample",
        "target name": "Gene",
        "ct": "Ct",
        "sample": "Sample",
        "gene": "Gene"
    })

    # Check required columns
    required = ["Sample", "Gene", "Ct"]

    missing_columns = [
        col for col in required
        if col not in raw_df.columns
    ]

    if missing_columns:
        st.error(
            "Missing columns: " +
            ", ".join(missing_columns)
        )
        st.stop()

    # Remove rows without sample or gene names
    raw_df = raw_df.dropna(
        subset=["Sample", "Gene"]
    ).copy()

    # Clean names
    raw_df["Sample"] = (
        raw_df["Sample"].astype(str).str.strip()
    )

    raw_df["Gene"] = (
        raw_df["Gene"].astype(str).str.strip()
    )

    # Convert Ct values to numbers
    # Undetermined becomes NaN (missing)
    raw_df["Ct"] = pd.to_numeric(
        raw_df["Ct"],
        errors="coerce"
    )

    if raw_df.empty:
        st.error("No valid sample rows were found.")
        st.stop()


    # ==========================================
    # 3. CALCULATE AVERAGE CT
    # ==========================================

    average_df = (
        raw_df.groupby(
            ["Sample", "Gene"],
            as_index=False
        )
        .agg(
            Average_Ct=("Ct", "mean"),
            Valid_Replicates=("Ct", "count"),
            Total_Replicates=("Ct", "size"),
            Ct_SD=("Ct", "std")
        )
    )

    average_df = average_df.rename(columns={
        "Average_Ct": "Average Ct",
        "Valid_Replicates": "Valid Replicates",
        "Total_Replicates": "Total Replicates",
        "Ct_SD": "Ct SD"
    })


    # ==========================================
    # 4. ANALYSIS SETTINGS
    # ==========================================

    st.header("Analysis Settings")

    # Let users choose a reference gene
    genes = sorted(raw_df["Gene"].unique())

    default_gene = (
        genes.index("GAPDH")
        if "GAPDH" in genes else 0
    )

    reference_gene = st.selectbox(
        "Select Reference Gene",
        options=genes,
        index=default_gene
    )

    # Let users choose normalization samples
    samples = sorted(
        sample for sample in raw_df["Sample"].unique()
        if sample.upper() != "NTC"
    )

    wt_defaults = [
        sample for sample in samples
        if sample in ["WT 1", "WT 2", "WT 3"]
    ]

    normalization_samples = st.multiselect(
        "Select Samples for Normalization",
        options=samples,
        default=wt_defaults
    )

    if not normalization_samples:
        st.warning(
            "Please select at least one "
            "normalization sample."
        )
        st.stop()


    # ==========================================
    # 5. PREPARE DATA FOR ANALYSIS
    # ==========================================

    # Exclude NTC from expression analysis
    analysis_df = average_df[
        average_df["Sample"].str.upper() != "NTC"
    ].copy()

    # Reference gene table
    reference_df = analysis_df[
        analysis_df["Gene"] == reference_gene
    ][
        ["Sample", "Gene", "Average Ct"]
    ].copy()

    reference_df = reference_df.rename(columns={
        "Gene": "Reference Gene",
        "Average Ct": "Reference Ct"
    })

    # Target gene table
    target_df = analysis_df[
        analysis_df["Gene"] != reference_gene
    ].copy()

    target_df = target_df.rename(columns={
        "Gene": "Target Gene",
        "Average Ct": "Target Ct"
    })

    if target_df.empty:
        st.error(
            "No target genes were found. "
            "Your file needs at least one gene "
            "besides the reference gene."
        )
        st.stop()

    # Match each sample to its reference Ct
    results = pd.merge(
        target_df,
        reference_df,
        on="Sample",
        how="left",
        validate="many_to_one"
    )


    # ==========================================
    # 6. CALCULATE DELTA CT
    # ==========================================

    results["Delta Ct"] = (
        results["Target Ct"] -
        results["Reference Ct"]
    )


    # ==========================================
    # 7. CALCULATE RELATIVE EXPRESSION
    # ==========================================

    results["Relative Expression"] = (
        2 ** (-results["Delta Ct"])
    )


    # ==========================================
    # 8. NORMALIZE TO SELECTED SAMPLES
    # ==========================================

    control_df = results[
        results["Sample"].isin(
            normalization_samples
        )
    ]

    control_average = (
        control_df.groupby(
            "Target Gene",
            as_index=False
        )["Relative Expression"]
        .mean()
    )

    control_average = control_average.rename(
        columns={
            "Relative Expression": "Control Average"
        }
    )

    normalized_df = pd.merge(
        results,
        control_average,
        on="Target Gene",
        how="left",
        validate="many_to_one"
    )

    normalized_df["Normalized Expression"] = (
        normalized_df["Relative Expression"] /
        normalized_df["Control Average"]
    )

    # Prevent invalid normalization results
    invalid_control = (
        normalized_df["Control Average"].isna() |
        (normalized_df["Control Average"] <= 0)
    )

    normalized_df.loc[
        invalid_control,
        "Normalized Expression"
    ] = float("nan")


    # ==========================================
    # 9. CREATE CLEAN RESULTS SUMMARY
    # ==========================================

    summary_df = normalized_df[
        [
            "Sample",
            "Target Gene",
            "Delta Ct",
            "Relative Expression",
            "Normalized Expression"
        ]
    ].copy()

    summary_df = summary_df.rename(columns={
        "Target Gene": "Gene",
        "Normalized Expression": "Normalized Fold"
    })


    # ==========================================
    # 10. DISPLAY WEBSITE TABS
    # ==========================================

    tab1, tab2, tab3 = st.tabs([
        "Raw Data",
        "Average Ct & QC",
        "Results"
    ])

    # TAB 1: RAW DATA
    with tab1:

        st.subheader("Raw qPCR Data")

        st.write(
            f"Total wells imported: {len(raw_df)}"
        )

        st.dataframe(
            raw_df,
            hide_index=True,
            use_container_width=True
        )


    # TAB 2: AVERAGE CT
    with tab2:

        st.subheader("Average Ct Results")

        st.dataframe(
            average_df.round(4),
            hide_index=True,
            use_container_width=True
        )

        missing_ct = raw_df["Ct"].isna().sum()

        st.write(
            f"Undetermined or missing Ct values: "
            f"{missing_ct}"
        )

        # Flag variable technical replicates
        # 0.5 Ct is an adjustable review threshold
        qc_flags = average_df[
            (average_df["Ct SD"] > 0.5) |
            (
                average_df["Valid Replicates"] <
                average_df["Total Replicates"]
            )
        ]

        if not qc_flags.empty:

            st.warning(
                "Some measurements need QC review."
            )

            st.dataframe(
                qc_flags.round(4),
                hide_index=True,
                use_container_width=True
            )

        else:
            st.success(
                "No replicate variability or "
                "missing-value flags detected."
            )

        st.caption(
            "QC flags are for review only. "
            "No Ct values are automatically removed. "
            "A passed check does not validate an assay."
        )


    # TAB 3: RESULTS
    with tab3:

        st.subheader("Gene Expression Results")

        st.write(
            f"Reference Gene: {reference_gene}"
        )

        st.write(
            "Normalization Samples: " +
            ", ".join(normalization_samples)
        )

        st.caption(
            "Normalization method: each sample's "
            "2^(-Delta Ct) divided by the arithmetic "
            "mean 2^(-Delta Ct) of the selected "
            "reference samples, separately for each gene."
        )

        st.dataframe(
            summary_df.round(6),
            hide_index=True,
            use_container_width=True
        )

        with st.expander(
            "View Complete Calculations"
        ):

            st.dataframe(
                normalized_df.round(6),
                hide_index=True,
                use_container_width=True
            )

        if invalid_control.any():
            st.warning(
                "Some genes have no valid "
                "normalization average. "
                "Review your selected reference "
                "samples and missing Ct values."
            )


        # ======================================
        # 11. EXPORT EXCEL WORKBOOK
        # ======================================

        buffer = BytesIO()

        with pd.ExcelWriter(
            buffer,
            engine="openpyxl"
        ) as writer:

            # Sheet 1
            raw_df.to_excel(
                writer,
                sheet_name="Raw Ct Data",
                index=False
            )

            # Sheet 2
            average_df.to_excel(
                writer,
                sheet_name="Average Ct",
                index=False
            )

            # Sheet 3
            normalized_df.to_excel(
                writer,
                sheet_name="Full Analysis",
                index=False
            )

            # Sheet 4
            summary_df.to_excel(
                writer,
                sheet_name="Results Summary",
                index=False
            )

            # Format all four sheets
            for sheet in writer.sheets.values():

                # Freeze headers
                sheet.freeze_panes = "A2"

                # Add filters
                sheet.auto_filter.ref = sheet.dimensions

                # Format headers
                for cell in sheet[1]:

                    cell.font = Font(
                        bold=True,
                        color="FFFFFF"
                    )

                    cell.fill = PatternFill(
                        fill_type="solid",
                        fgColor="24476B"
                    )

                # Center all cells
                for row in sheet.iter_rows():

                    for cell in row:

                        cell.alignment = Alignment(
                            horizontal="center",
                            vertical="center"
                        )

                        if isinstance(
                            cell.value, float
                        ):
                            cell.number_format = "0.000000"

                # Adjust column widths
                for column in sheet.columns:

                    letter = get_column_letter(
                        column[0].column
                    )

                    max_length = max(
                        len(str(cell.value))
                        if cell.value is not None else 0
                        for cell in column
                    )

                    sheet.column_dimensions[
                        letter
                    ].width = min(
                        max(max_length + 3, 16),
                        32
                    )

        # Excel download button
        st.download_button(
            label="Download Excel Results",
            data=buffer.getvalue(),
            file_name="PlateLab_Results.xlsx",
            mime=(
                "application/vnd.openxmlformats-officedocument."
                "spreadsheetml.sheet"
            )
        )
