import streamlit as st

import_and_prepare = st.Page(
    page="page_1.py",
    title="Import and parametrisation",
)

calibration = st.Page(
    page="page_2.py",
    title="Calibration",
)

pg = st.navigation(
    [
        import_and_prepare,
        calibration,
    ]
)

pg.run()
