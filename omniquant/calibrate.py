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
import math
from pathlib import Path
from collections import namedtuple

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

    def __init__(self, name, x, y, case, degree=2, weight=None):

        # Private
        self._polynome = None
        self._polynome_plot = None

        # Public
        self.name = name
        self.x = x
        self.y = y
        self.case = case
        self.degree = degree
        self.weight = weight
        self.excluded_values = {
            "x": [],
            "y": []
        }
        self.polynome_details = None

        match weight:
            case "1/X":
                self.weight = 1 / math.sqrt(self.x)
            case "1/X²":
                self.weight = 1 / self.x
            case _:
                self.weight = None

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
        """
        Function to reset the polynomials. Should be called when any method modifies the calibration data (exclusion of
        some data points for example)
        """
        self._scaled_polynome = None
        self._polynome_plot = None

    def drop(self, axis, value):
        """
        Drop a value and associated value from axes x and y using index. The index and the value are saved to the
        excluded values.
        :param axis: axis x or y on which to search for value to remove. The removal of a value from an axis removes
                     its sister value from the other axis.
        :param value: value to remove.
        """

        # TODO: Use masks on arrays to keep track of points to use

        if axis not in ["x", "y"]:
            raise ValueError(f"Axis term can only be x or y. Detected term: {axis}")
        # Get index where value to drop is situated
        idx = np.where(getattr(self, axis) == value)[0]
        # keep record of deleted values and associated index
        self.excluded_values["x"].append((idx, self.x[idx]))
        self.excluded_values["y"].append((idx, self.y[idx]))
        self.x, self.y = np.delete(self.x, idx), np.delete(self.y, idx)
        self._reset()

    @property
    def equation(self):
        """
        Return's the equation by converting the [-1, 1] scaled polynomial to the original scale.
        :return: Converted polynomial to original scale
        """
        # Get the converted equation coefficients
        coefficients = list(self.scaled_polynome.convert().coef)
        # Reverse to get the coefficients in the standard order
        coefficients.reverse()
        return np.poly1d(coefficients)

    @property
    def scaled_polynome(self):
        """
        Get the scaled polynome from the Polynomial Fit class of numpy
        :return: Scaled ([-1, 1]) polynome.
        """

        if self._polynome:
            return self._polynome
        # Create polynome w/ calibration data (x=concentrations, y=signals/ratios)
        self._polynome, polynome_details = np.polynomial.Polynomial.fit(
            x=self.x,
            y=self.y,
            deg=self.degree,
            full=True,
            w=self.weight  # TODO: test ponderation compared to commercial software
        )
        self.polynome_details = {
            "SSR": polynome_details[0],
            "Rank": polynome_details[1],
            "SV": polynome_details[2],
            "Rcond": polynome_details[3]
        }
        return self._polynome

    @property
    def polynome_plot(self):
        """
        Get the updated polynomial plot from the Plotly interface.
        :return: Plotly figure
        """

        if self._polynome_plot:
            return self._polynome_plot

        x_name = "Metabolite concentration"
        match self.case:
            case 2:
                y_name = "Signal"
            case 4:
                y_name = "Signal/IS"
            # Wildcard case to catch other cases as an error
            case _:
                raise ValueError("Couldn't coerce name for y axis. Please make sure you have calibration curve data")

        # Generate the graph object containing the calibration points
        points = go.Scatter(
            x=self.x,
            y=self.y,
            mode='markers',
            name='data'
        )
        # Generate fitted curve to calibration points
        poly_x, poly_y = self.scaled_polynome.linspace()
        trend = go.Scatter(
            x=poly_x,
            y=poly_y,
            mode='lines'
        )
        # Draw out figure
        data = [points, trend]
        self._polynome_plot = go.Figure(data=data)
        self._polynome_plot.update_layout({
            "title": f'{self.name}',
            "showlegend": True,
            "xaxis_title": x_name,
            "yaxis_title": y_name
        }
        )
        return self._polynome_plot

    @property
    def residuals(self):
        if self._residuals:
            return self._residuals
        return (self.equation(self.x) - self.y) / self.y

    @property
    def limits(self):
        """
        Lower and upper limits for the calibration curve
        :return: namedtuple containing limits
        """

        Limits = namedtuple("Limits", "lower upper")
        return Limits(self.y.min(), self.y.max())


class Quantifier:
    CASES = [1, 2, 3, 4, 5]

    def __init__(
            self,
            metabolite_name: str,
            cal_data: pd.DataFrame
    ):

        # Private
        self._calibrator = None
        self._case = None

        # Public
        self.metabolite_name = metabolite_name
        self.cal_data = cal_data
        self.is_int_std = False
        self.is_int_std_conc_known = False
        self.is_cal_points = False
        self.quantify = None
        self.response_factor = 2
        self.case = self.select_case()

    def __repr__(self):

        return (f"Metabolite Name: {self.metabolite_name}\n"
                f"Calibration Data: {self.cal_data}\n"
                f"Is an internal standard present: {'Yes' if self.is_int_std else 'No'}\n"
                f"Is the internal standard concentration known: {'Yes' if self.is_int_std_conc_known else 'No'}\n"
                f"Are there multiple calibration points: {'Yes' if self.is_cal_points else 'No'}")

    @property
    def case(self):
        return self._case

    @case.setter
    def case(self, value):
        if value not in self.CASES:
            raise ValueError(f"Case possible choice are: {self.CASES}")
        self._case = value
        self._set_quantifier_func()

    def select_case(self):
        """
        Return the use case we are in

        :return: Use case
        """

        # TODO: Function needs refactoring (too nested for nothing)
        case = 1
        if not self.cal_data.empty:
            if self.cal_data["Cal_Signal"].any() and len(self.cal_data["Cal_Signal"] > 1):
                self.is_cal_points = True
                case = 2
            if self.cal_data["IS_signal"].any():
                self.is_int_std = True
                if not self.cal_data["IS_Concentration"].any():
                    case = 3
                    if self.is_cal_points:
                        case = 4
                else:
                    self.is_int_std_conc_known = True
                    case = 5
        return case

    def _set_quantifier_func(self):
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

    @property
    def calibrator(self):
        """
        Initialize or get the calibrator. Only used in cases 2 & 4.
        :return: Initialized calibrator object
        """

        if self._calibrator:
            return self._calibrator
        if self.case not in [2, 4]:
            raise CaseError(f"Use case is initialized as {self.case}. Only cases 2 & 4 get a calibrator")

        try:
            self._calibrator = Calibrator(
                name=self.metabolite_name,
                x=self.cal_data["Cal_Concentration"].to_numpy(),
                y=self.cal_data["Cal_Signal"].to_numpy() if self.is_int_std is False
                else np.divide(self.cal_data["Cal_Signal"].to_numpy(), self.cal_data["IS_signal"]),
                case=self.case
            )
        except Exception:
            raise CalibrationError("There was an error while initializing the carlibrator")

        return self._calibrator

    @staticmethod
    def apply_response_factor(func):
        """
        Wrapper around quantification function to apply the response factor
        :param func: function to wrap (self.quantify)
        :return: wrapped function
        """

        def wrapper(self, *args, **kwargs):
            return func(self, *args, **kwargs) * self.response_factor

        return wrapper

    @staticmethod
    def _quantify_no_int_std_no_curve(signal):
        """
        Relative quantification of our signal.
        :param signal: Metabolite sample signal/area
        :return: Metabolite sample signal/area
        """
        if not isinstance(signal, float) or not isinstance(signal, int):
            raise ValueError(f"Signal should be a number. Detected type: {type(signal)}:")
        return signal

    def _quantify_no_int_std_with_curve(self, signal):
        return (self.calibrator.equation - signal).roots[1]

    @staticmethod
    def _quantify_int_std_no_conc_no_curve(signal, std_signal):
        return signal / std_signal

    def _quantify_int_std_no_conc_with_curve(self, signal, std_signal):
        to_quant = signal / std_signal
        if to_quant > self.calibrator.limits.upper:
            raise CalibrationError("Under calibration range")
        if to_quant < self.calibrator.limits.lower:
            raise CalibrationError("Over calibration range")
        return (self.calibrator.equation - to_quant).roots[1]

    @apply_response_factor
    def _quantify_int_std_with_conc(self, signal):
        return signal * self.cal_data["IS_Concentration"]


class CaseError(Exception):
    pass


class CalibrationError(Exception):
    pass


class QuantificationError(Exception):
    pass


if __name__ == "__main__":
    pd.options.display.max_columns = None

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
    test_case_5 = (
        "Acetate",
        pd.DataFrame.from_dict(
            {
                "Cal_name": "",
                "Cal_Signal": "",
                "Cal_Concentration": "",
                "IS_signal": np.array([0.703773953]),
                "IS_Concentration": np.array([4.2])
            }
        )
    )

    test_data = [
        test_case_1,
        test_case_2,
        test_case_3,
        test_case_4,
        test_case_5
    ]

    quant = Quantifier(test_case_4[0], test_case_4[1])
    print(quant.calibrator.scaled_polynome)
    # print(list(quant.calibrator.scaled_polynome.convert().coef).reverse())
    quant.calibrator.polynome_plot.show()
    print(quant.calibrator.equation)
    print(quant.quantify(1163372160, 483495840))
    print((quant.calibrator.equation(quant.calibrator.x) - quant.calibrator.y) / quant.calibrator.y)
    print(quant.calibrator.limits)
    # print(f"a = {a}")
    # a.reverse()
    # print(f"a2 = {a}")
    # b = np.poly1d(a)
    # print(f"b = {b}")
    # conc = (b - 70279200).roots
    # print(conc)
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
