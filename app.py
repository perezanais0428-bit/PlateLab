
# PLATELAB - Our First Website

import streamlit as st
import pandas as pd

# Website title
st.title("🧬 PlateLab")

st.write("qPCR Experiment Design & Analysis")

st.write(
    "Upload your qPCR results to automatically "
    "calculate average Ct values."
)

# Create a file upload button
uploaded_file = st.file_uploader(
    "Upload qPCR CSV file",
    type=["csv"]
)

# Only analyze data after a file is uploaded
if uploaded_file is not None:

    # Read uploaded data
    raw_df = pd.read_csv(uploaded_file)

    # Standardize column names
    raw_df.columns = raw_df.columns.str.strip().str.lower()

    raw_df = raw_df.rename(columns={
        "sample name": "Sample",
        "target name": "Gene",
        "ct": "Ct"
    })

    # Check the required columns
    required = ["Sample", "Gene", "Ct"]

    if not all(col in raw_df.columns for col in required):
        st.error("Missing required CSV columns.")
        st.stop()

    # Convert Ct values to numbers
    raw_df["Ct"] = pd.to_numeric(
        raw_df["Ct"], errors="coerce"
    )

    st.subheader("Raw qPCR Data")
    st.dataframe(raw_df, hide_index=True)

    # Calculate average Ct
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

    st.subheader("Average Ct Results")

    st.dataframe(
        average_df.round(4),
        hide_index=True
    )

    st.success("Your qPCR data was imported!")
