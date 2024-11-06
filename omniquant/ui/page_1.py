import tkinter as tk
from pathlib import Path
from tkinter import filedialog

import pandas as pd
import streamlit as st
from sess_i.base.main import SessI

from omniquant import __version__
from omniquant.io.dataset import CalibrationData, SampleData
from omniquant.numeric.quantify import Quantifier


#############
# FUNCTIONS #
#############

@st.cache_data
def read_sample_data(path: str) -> SampleData:
    if path:
        data = SampleData.from_file(path)
        return data


@st.cache_data
def read_cal_data(path: str) -> CalibrationData:
    if path:
        data = CalibrationData.from_file(path)
        return data


# noinspection PyStatementEffect
def set_state(state: int):
    # noinspection PyStatementEffect
    st.session_state[f"page_1_state"] = state


def select_file() -> str:
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


def generate_quantifiers(dataset: SampleData) \
        -> tuple[dict, dict]:
    """
    Generate quantifiers for the dataset.

    This function takes a SampleData object and generates two lists of
    Quantifier objects. One dict contains quantifiers with calibration data,
    and the other dict contains quantifiers without calibration data.

    :param dataset: SampleData object containing the dataset
    :return: A tuple containing two dicts of quantifier objects
    """
    no_cal_quants = {}
    cal_quants = {}
    for analyte in dataset.analytes:
        quant_obj = Quantifier(
            metabolite_name=analyte,
            cal_data=dataset.cal_points[
                dataset.cal_points.Analyte == analyte
                ],
        )
        if quant_obj.case in [2, 4]:
            cal_quants[analyte] = quant_obj
        else:
            no_cal_quants[analyte] = quant_obj
    return cal_quants, no_cal_quants


########
# MAIN #
########


st.set_page_config(
    page_title="OmniQuant",
    page_icon=":bar_chart:",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title(f"OmniQuant v{__version__}")

# Initialize session
session = SessI(session_state=st.session_state)

# Build sidebar
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

st.sidebar.divider()
st.sidebar.title("Options")
reset_process = st.sidebar.button(
    label="Reset current process"
)
st.sidebar.divider()

# Input data
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

submit_data = st.button(
    label="Submit Data",
    on_click=set_state,
    args=[0]
)

if st.session_state.get("page_1_state") == 0:
    # Display data and allow for selection of metabolites to drop
    with st.expander(label="metabolite_selector", expanded=False):
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
        st.button(label="Submit", on_click=set_state, args=[1])

    # Generate quantifiers
    if st.session_state.get("page_1_state") == 1:
        st.write(sample_data.cal_points)
        with_cal = session.get_object("with_cal")
        no_cal = session.get_object("no_cal")
        if not with_cal and not no_cal:
            with_cal, no_cal = generate_quantifiers(sample_data)
            session.register_object(with_cal, key="with_cal")
            session.register_object(no_cal, key="no_cal")

        # parametrise quantification
        with st.form(key="quantification_form"):
            def save_form_values():
                st.session_state.weight = st.session_state.weight_select
                st.session_state.regression = st.session_state.regression_select

            if with_cal:
                weight_index = {
                    "None": 0,
                    "1/X": 1,
                    "1/X²": 2
                }
                weight = st.selectbox(
                    label="Select weight",
                    options=["None", "1/X", "1/X²"],
                    key="weight_select",
                    index=weight_index[st.session_state.weight]
                    if "weight" in st.session_state else 0
                )
                regress_index = {
                    "None": 0,
                    "Linear": 1,
                    "Quadratic": 2
                }
                regression = st.selectbox(
                    label="Select regression",
                    options=["None", "Linear", "Quadratic"],
                    key="regression_select",
                    index=regress_index[st.session_state.regression]
                    if "regression" in st.session_state else 0
                )

            st.form_submit_button(
                label="Submit Parameters", on_click=save_form_values
            )
