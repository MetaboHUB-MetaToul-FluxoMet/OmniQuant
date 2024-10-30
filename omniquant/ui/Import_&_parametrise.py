import tkinter as tk
from tkinter import filedialog
from pathlib import Path

import pandas as pd
import streamlit as st

from omniquant import __version__
from omniquant.numeric.quantify import Quantifier
from omniquant.io.dataset import CalibrationData, SampleData
from sess_i.base.main import SessI


@st.cache_data
def read_sample_data(path):
    if path:
        data = SampleData.from_file(path)
        return data


@st.cache_data
def read_cal_data(path):
    if path:
        data = CalibrationData.from_file(path)
        return data


def select_file():
    # Set up tkinter for directory chooser
    root = tk.Tk()
    root.withdraw()

    # Make folder picker dialog appear on top of other windows
    root.wm_attributes('-topmost', 1)

    return str(Path(
        filedialog.askopenfilename(
            master=root,
            title="Select the data file",
            filetypes=[("Data files", "*.tsv *.txt")]
        )
    ))


st.set_page_config(
    page_title="OmniQuant",
    page_icon=":bar_chart:",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title(f"OmniQuant v{__version__}")

# Initialize session and process
session = SessI(session_state=st.session_state)

st.sidebar.title("Ressources")
st.sidebar.link_button(
    label="Documentation",
    url="https://multinmrfit.readthedocs.io/en/latest/",
    disabled=True
)
st.sidebar.link_button(
    label="GitHub project",
    url="https://github.com/llegregam/OmniQuant"
)

# Build sidebar
st.sidebar.divider()
st.sidebar.title("Options")
reset_process = st.sidebar.button(
    label="Reset current process"
)
st.sidebar.divider()


with st.container():
    st.subheader("Import Data")
    btn_cols = st.columns(2)

    # Initialize buttons for selecting data
    with btn_cols[0]:
        sample_select_btn = st.button("Select Sample Data")
    with btn_cols[1]:
        cal_select_btn = st.button("Select Calibration Data")

# Handle paths
if sample_select_btn:
    sample_data_path = select_file()
    sample_data = read_sample_data(sample_data_path)
    session.register_object(sample_data, key="sample_data")

if cal_select_btn:
    cal_data_path = select_file()
    cal_data = read_cal_data(cal_data_path)
    session.register_object(cal_data, key="cal_data")

sample_data = session.get_object("sample_data")
cal_data = session.get_object("cal_data")

if sample_data:
    # Display data
    with st.expander(label="Metabolite Selector", expanded=True):
        to_drop = st.multiselect(
            label="Select Metabolites to drop",
            options=sample_data.analytes if sample_data else []
        )
        sample_data.data = sample_data()[~sample_data().Analyte.isin(to_drop)]
        st.dataframe(sample_data())
        if cal_data:
            cal_data.data = cal_data()[~cal_data().Analyte.isin(to_drop)]
            st.dataframe(cal_data())

        sample_data.data = pd.merge(
            left=sample_data(),
            right=cal_data(),
            how="left",
        )
        st.write(sample_data.cal_points)

    # parametrise quantification
    with st.expander(label="Quantification Parameters", expanded=True):
        for analyte in sample_data.analytes:
            quantifier = Quantifier(
                metabolite_name=analyte,
                cal_data=sample_data.cal_points[
                    sample_data.cal_points.Analyte == analyte
                ],
            )
            session.register_object(quantifier, key=f"{analyte}_quantifier")
            st.write(quantifier)
            st.write(f"Case = {quantifier.case}")
            if quantifier.case in [2, 4]:
                st.write(quantifier.calibrator.polynome_plot)
            st.write(quantifier.cal_data)

