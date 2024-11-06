import streamlit as st
from sess_i.base.main import SessI

session = SessI(session_state=st.session_state)

# st.write(
#     st.session_state
# )
if not session.get_object("with_cal"):
    st.header("❌No calibration data available❌")
else:
    st.header("✅Calibration data available✅")
    with_cal = session.get_object("with_cal")
    for analyte in with_cal:
        with (st.expander(f"Quantify {analyte}")):
            quantifier = with_cal[analyte]
            st.write(quantifier)
            st.write(quantifier.calibrator.polynome_plot)

