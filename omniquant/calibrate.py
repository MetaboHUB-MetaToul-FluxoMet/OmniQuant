"""
Quantifier module for OmniQuant.

The Quantifier's goal is to build select which use case is present following the decision tree diagram. It dynamically
selects which quantify method should be called, and decides if a calibration curve is needed. It also computes the
polynomials, uses it to fit the data and make predictions for quantification. It can handle multiple types of internal
standards (C12 and C13 for example).

Steps:
* From the calibration data dataframe, parse information and select use case between the 5 available cases
* Set's up the quantify function
* Builds the calibration and the calibration curves if needed
* Runs quantification
* Returns the calculated data

"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objs as go


def read_data(path: str) -> pd.DataFrame:
    if not isinstance(path, str):
        raise TypeError(f"Path to data must be a string. Detected type: {type(path)}")
    data_path = Path(path)
    if not data_path.exists():
        raise ValueError(f"The given data path does not exist. Path: {str(data_path)}")
    df = pd.read_csv(data_path, sep=",")
    return df


class Calibrator:

    def __init__(self, name, x, y, case, degree=2):

        self.name = name
        self.x = x
        self.y = y
        self.case = case
        self.degree = degree
        self._polynome = None
        self._polynome_plot = None
        self.excluded_values = {
            "x" : [],
            "y" : []
        }

    def __setattr__(self, key, value):

        # We handle x and y setting here because the same verification procedure is applied to both.
        if key in ["x", "y"]:
            if type(value) != np.ndarray:
                try:
                    value = np.array(value)
                except Exception:
                    raise TypeError(f"There was an error while converting the array to numpy ndarray. Array: {value}")
            if len(value) <= 1:
                raise ValueError(f"To build calibration, x values must be more than 1. Number of values: {len(value)}")
            if not value.all() > 0:
                raise ValueError(f"All the values of the calibration data arrays must be positive")

        if (key == "case") and (value not in [2, 4]):
            raise ValueError(f"Use case can only be 2 or 4. Detected case: {key}")

        self.__dict__[key] = value

    def _reset(self):

        self.scaled_polynome = None
        self._polynome_plot = None
        self._equation = None

    def drop_value(self, axis, value):

        if axis not in ["x", "y"]:
            raise ValueError(f"Axis term can only be x or y. Detected term: {axis}")

        idx = np.where(getattr(self, axis) == value)[0]
        self.excluded_values["x"].append((idx, self.x[idx]))
        self.excluded_values["y"].append((idx, self.y[idx]))
        self.x = np.delete(self.x, idx)
        self.y = np.delete(self.y, idx)
        self._reset()

    @property
    def equation(self):

        return self.scaled_polynome.convert()

    @property
    def scaled_polynome(self):

        if self._polynome:
            return self._polynome
        self._polynome = np.polynomial.Polynomial.fit(
            x=self.x,
            y=self.y,
            deg=self.degree
        )
        return self._polynome

    @property
    def polynome_plot(self):

        if self._polynome_plot:
            return self._polynome_plot

        x_name = "Metabolite concentration"
        y_name = None
        match self.case:
            case 2:
                y_name = "Signal"
            case 4:
                y_name = "Signal/IS"
        if not y_name:
            raise ValueError("Couldn't coerce name for y axis. Please make sure you have calibration curve data")

        points = go.Scatter(
            x=self.x,
            y=self.y,
            mode='markers',
            name='data'
        )
        poly_x, poly_y = self.scaled_polynome.linspace()
        trend = go.Scatter(
            x=poly_x,
            y=poly_y,
            mode='lines'
        )
        layout = go.Layout(
            title=f'{self.name}',
            showlegend=True,
            xaxis_title=x_name,
            yaxis_title=y_name
        )
        data = [points, trend]
        self._polynome_plot = go.Figure(data=data, layout=layout)
        return self._polynome_plot


class Quantifier:

    def __init__(
            self,
            metabolite_name: str,
            cal_data: pd.DataFrame
    ):
        self.metabolite_name = metabolite_name
        self.cal_data = cal_data
        self.is_int_std = False
        self.is_int_std_conc_known = False
        self.is_cal_points = False
        self.quantify = None
        self.case = 1
        if not self.cal_data.empty:
            if self.cal_data["Cal_Signal"].any() and len(self.cal_data["Cal_Signal"] > 1):
                self.is_cal_points = True
                self.case = 2
            if self.cal_data["IS_signal"].any():
                self.is_int_std = True
                if not self.cal_data["IS_Concentration"].any():
                    self.case = 3
                    if self.is_cal_points:
                        self.case = 4
                else:
                    self.is_int_std_conc_known = True
                    self.case = 5
        self.calibrator = Calibrator(
            name=self.metabolite_name,
            x=cal_data["Cal_Concentration"].to_numpy(),
            y=self.cal_data["Cal_Signal"].to_numpy() if self.is_int_std is False
            else np.divide(self.cal_data["Cal_Signal"].to_numpy(), self.cal_data["IS_signal"]),
            case=self.case
        )
        self.set_quantifier()

    def __repr__(self):

        return (f"Metabolite Name: {self.metabolite_name}\n"
                f"Calibration Data: {self.cal_data}\n"
                f"Is an internal standard present: {'Yes' if self.is_int_std else 'No'}\n"
                f"Is the internal standard concentration known: {'Yes' if self.is_int_std_conc_known else 'No'}\n"
                f"Are there multiple calibration points: {'Yes' if self.is_cal_points else 'No'}")

    def set_quantifier(self):
        match self.case:
            case 1:
                self.quantify = self._quantify_no_int_std_no_curve
            case 2:
                self.quantify = self._quantify_no_int_std_with_curve
            case 3:
                self.quantify = self._quantify_int_std_no_conc_no_curve
            case 4:
                self.quantify = self._quantify_int_std_no_conc_with_curve
            case 5:
                self.quantify = self._quantify_int_std_with_conc

    def _quantify_no_int_std_no_curve(self, signal):
        """
        Relative quantification of our signal.
        :param signal: Metabolite sample signal/area
        :return: Metabolite sample signal/area
        """
        if not isinstance(signal, float) or not isinstance(signal, int):
            raise ValueError(f"Signal should be a number. Detected type: {type(signal)}:")
        return signal

    def _quantify_no_int_std_with_curve(self):
        return "Case 2"

    def _quantify_int_std_no_conc_no_curve(self):
        return "Case 3"

    def _quantify_int_std_no_conc_with_curve(self):
        return "Case 4"

    def _quantify_int_std_with_conc(self):
        return "Case 5"


if __name__ == "__main__":
    pd.options.display.max_columns = None
    # COLS = [
    #     "File Name",
    #     "Sample Type",
    #     "Molecule",
    #     "Total Area",
    #     "Internal Standard Concentration"
    # ]
    # data = read_data(r"C:\Users\legregam\PycharmProjects\OmniQuant\omniquant\test-data\test_data.csv")

    # Case 1: No IS, no cal points
    test_case_1 = (
        "ATP",
        pd.DataFrame.from_dict(
            {
                "Cal_name": [],
                "Cal_Signal": None,
                "Cal_Concentration": None,
                "IS_signal": None,
                "IS_Concentration": None
            }
        )
    )
    # Case 2: No IS, calibration points
    test_case_2 = (
        "dATP",
        pd.DataFrame.from_dict(
            {
                "Cal_name": ["SMS_01", "SMS_02", "SMS_04", "SMS_08", "SMS_16", "SMS_31"],
                "Cal_Signal": np.array([89045904, 83217648, 70279200, 65918260, 37007596, 17540828]),
                "Cal_Concentration": np.array([2.618, 1.309, 0.655, 0.327, 0.164, 0.082]),
                "IS_signal": None,
                "IS_Concentration": None
            }
        )
    )
    # Case 3: IS present, No [IS], no cal points
    test_case_3 = (
        "ATP",
        pd.DataFrame.from_dict(
            {
                "Cal_name": ["SMS_01", "SMS_02", "SMS_04", "SMS_08", "SMS_16", "SMS_31"],
                "Cal_Signal": None,
                "Cal_Concentration": None,
                "IS_signal": np.array([483495840, 436619616, 412101248, 381654080, 345984544, 316497536]),
                "IS_Concentration": None
            }
        )
    )
    # Case 4: IS present, no [IS], calibration points
    test_case_4 = (
        "ADP",
        pd.DataFrame.from_dict(
            {
                "Cal_name": ["SMS_01", "SMS_02", "SMS_04", "SMS_08", "SMS_16", "SMS_31"],
                "Cal_Signal": np.array([1163372160, 506132736, 234740656, 107778792, 48474812, 22390422]),
                "Cal_Concentration": np.array([2.618, 1.309, 0.655, 0.327, 0.164, 0.082]),
                "IS_signal": np.array([483495840, 436619616, 412101248, 381654080, 345984544, 316497536]),
                "IS_Concentration": None
            }
        )
    )
    # Case 5: IS present, [IS] known (NMR with TSP for example)

    test_data = [
        test_case_1,
        test_case_2,
        test_case_3,
        test_case_4
    ]

    quant = Quantifier(test_case_2[0], test_case_2[1])
    print(quant.calibrator.scaled_polynome)
    print(quant.calibrator.equation)
    quant.calibrator.polynome_plot.show()

    # print(data.head(10))
    # for metabolite in data.Molecule.unique():
    #     met_data = data.loc[data["Molecule"] == metabolite]
    #     met_data = met_data[COLS]
    #     quantifier = Quantifier(
    #         metabolite_name=metabolite,
    #         cal_data=met_data.loc[met_data["Sample Type"] == "Standard", "Total Area"].values,
    #         signal=met_data.loc[met_data["Sample Type"] == "Unknown", "Total Area"].values,
    #         int_std="IDMS",
    #         int_std_conc=None,
    #         is_int_std=True,
    #         is_int_std_conc_known=False,
    #         is_cal_points=True
    #     )
    #     print(quantifier)
