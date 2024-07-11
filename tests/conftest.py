import pytest
import pandas as pd
import numpy as np

from omniquant.calibrate import Quantifier

@pytest.fixture(params=[
    # Case 1: No IS, no cal points
    ("ATP", pd.DataFrame.from_dict({
        "Cal_name": [],
        "Cal_Signal": None,
        "Cal_Concentration": None,
        "IS_Signal": None,
        "IS_Concentration": None
    })),
    # Case 2: No IS, calibration points
    ("dATP", pd.DataFrame.from_dict({
        "Cal_name": ["SMS_01", "SMS_02", "SMS_04", "SMS_08", "SMS_16", "SMS_31"],
        "Cal_Signal": np.array([89045904, 83217648, 70279200, 65918260, 37007596, 17540828]),
        "Cal_Concentration": np.array([2.618, 1.309, 0.655, 0.327, 0.164, 0.082]),
        "IS_Signal": None,
        "IS_Concentration": None
    })),
    # Case 3: IS present, No [IS], no cal points
    ("ATP", pd.DataFrame.from_dict({
        "Cal_name": ["SMS_01", "SMS_02", "SMS_04", "SMS_08", "SMS_16", "SMS_31"],
        "Cal_Signal": None,
        "Cal_Concentration": None,
        "IS_Signal": np.array([483495840, 436619616, 412101248, 381654080, 345984544, 316497536]),
        "IS_Concentration": None
    })),
    # Case 4: IS present, no [IS], calibration points
    ("ADP", pd.DataFrame.from_dict({
        "Cal_name": ["SMS_01", "SMS_02", "SMS_04", "SMS_08", "SMS_16", "SMS_31"],
        "Cal_Signal": np.array([1163372160, 506132736, 234740656, 107778792, 48474812, 22390422]),
        "Cal_Concentration": np.array([2.618, 1.309, 0.655, 0.327, 0.164, 0.082]),
        "IS_Signal": np.array([483495840, 436619616, 412101248, 381654080, 345984544, 316497536]),
        "IS_Concentration": None
    })),
    # Case 5: IS present, [IS] known (NMR with TSP for example)
    ("Acetate", pd.DataFrame.from_dict({
        "Cal_name": "",
        "Cal_Signal": "",
        "Cal_Concentration": "",
        "IS_Signal": np.array([0.703773953]),
        "IS_Concentration": np.array([4.2])
    }))
])
def quantifier_instance(request):
    metabolite_name, cal_data = request.param
    return Quantifier(metabolite_name, cal_data, calib_weight="1/X²")

