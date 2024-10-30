"""
Quantifier module for OmniQuant.

The Quantifier's goal is to build select which use case is present following
the decision tree diagram. It dynamically selects which quantify method
should be called, and decides if a calibration curve is needed. It also
computes the polynomials, uses it to fit the data and make predictions for
quantification. It can handle multiple types of internal standards (C12 and
C13 for example).

Steps: * From the calibration data dataframe, parse information and select
use case between the 5 available cases * Set's up the quantify function *
Builds the calibration and the calibration curves if needed * Runs
quantification * Returns the calculated data

"""
from collections import namedtuple
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import plotly.graph_objs as go


def read_data(path: str) -> pd.DataFrame:
    """
    This function reads a CSV file from a given path and returns it as a
    pandas DataFrame.

    Parameters:
    path (str): The path to the CSV file.

    Returns:
    pd.DataFrame: The DataFrame created from the CSV file.
    """
    return pd.read_csv(Path(path), sep=",")


class Calibrator:

    def __init__(self, name, x, y, case, weight, degree=2):

        # Private
        self._polynome: Union[np.polynomial.Polynomial, None] = None
        self._polynome_plot: Union[go.Figure, None] = None
        self._residual_plot: Union[go.Figure, None] = None
        self._x: np.ndarray = x
        self._y: np.ndarray = y
        self._mask: Union[np.ndarray, np.ma.MaskedArray] = np.full(
            shape=self._x.shape,
            fill_value=False)

        # Public
        self.name: str = name
        self.case: int = case
        self.degree: int = degree
        self.excluded_values: dict = {
            "x": np.ma.masked_array(self._x, mask=~self._mask).compressed(),
            "y": np.ma.masked_array(self._y, mask=~self._mask).compressed()
        }
        self.polynome_details: Union[np.ndarray, None] = None

        # Weighting used for the calibration curve
        weight_dict: dict = {
            "1/X": np.full_like(self.x, 1) / self.x,
            "1/X²": np.full_like(self.x, 1) / self.x ** 2
        }
        self.weighted_x: np.ndarray = weight_dict.get(weight)

    def __setattr__(self, key: str, value: Union[np.ndarray, None]):

        # We handle x and y setting here because the same verification
        # procedure is applied to both.
        if key in ["x", "y"]:
            # Must be numpy array
            value = np.array(value) if not isinstance(value,
                                                      np.ndarray) else value
            if len(value) <= 1:
                raise ValueError(
                    f"To build calibration, x values must be more than 1. "
                    f"Number of values: {len(value)}")
            if not value.all() > 0:
                raise ValueError(
                    f"All the values of the calibration data arrays must be "
                    f"positive")

            self.__dict__[f"_{key}"] = value

        elif key == "case" and value not in [2, 4]:
            raise ValueError(
                f"Use case can only be 2 or 4. Detected case: {value}")

        else:
            self.__dict__[key] = value

    # x and y axes are handled using a common mask to ease data point
    # exclusion and inclusion
    @property
    def x(self) -> np.ndarray:
        return np.ma.masked_array(self._x, mask=self._mask).compressed()

    @property
    def y(self) -> np.ndarray:
        return np.ma.masked_array(self._y, mask=self._mask).compressed()

    def update_axes(self, indice: int):
        """
        Update the axes mask at given indice and make value opposite

        :param indice: indice at which the given value should be included/excluded
        """

        self._mask[indice] = not self._mask[indice]
        self._reset()

    def _reset(self):
        """
        Function to reset the polynomials. Should be called when any method modifies the calibration data (exclusion of
        some data points for example)
        """
        self._polynome = None
        self._scaled_polynome = None
        self._polynome_plot = None


    @property
    def equation(self) -> np.poly1d:
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
    def scaled_polynome(self) -> np.polynomial.Polynomial:
        """
        Get the scaled polynome from the Polynomial Fit class of numpy
        :return: Scaled ([-1, 1]) polynome.
        """

        if self._polynome:
            return self._polynome

        # Fit the polynome to the data
        self._polynome, polynome_details = np.polynomial.Polynomial.fit(
            x=self.x,
            y=self.y,
            deg=self.degree,
            full=True,
            w=self.weighted_x
        )
        self.polynome_details = {
            "SSR": polynome_details[0],
            "Rank": polynome_details[1],
            "SV": polynome_details[2],
            "Rcond": polynome_details[3]
        }
        return self._polynome

    @property
    def polynome_plot(self) -> go.Figure:
        """
        Get the updated polynomial plot from the Plotly interface.
        :return: Plotly figure
        """

        if self._polynome_plot:
            return self._polynome_plot

        x_name = "Metabolite concentration"
        y_name_dict = {2: "Signal", 4: "Signal/IS"}
        if self.case not in y_name_dict:
            raise ValueError(
                "Couldn't coerce name for y axis. Please make sure you have "
                "calibration curve data")
        y_name = y_name_dict[self.case]

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
    def residuals(self) -> np.ndarray:
        """
        Calculation of the residuals for the calibration curve. Observed
        values are self.y and predicted values are self.equation(self.x).
        The residuals are calculated as the difference between the observed
        and predicted values divided by the observed values to normalize the
        residuals.

        :return: vector containing residuals
        """
        return (self.equation(self.x) - self.y) / self.y

    @property
    def residual_plot(self) -> go.Figure:
        """
        Draw the residuals plot
        :return: plotly figure
        """

        if self._residual_plot:
            return self._residual_plot

        residuals = self.residuals
        residuals_plot = go.Figure()
        residuals_plot.add_trace(
            go.Scatter(x=self.x, y=residuals, mode='markers'))
        residuals_plot.update_layout(
            title=f"Residuals plot for {self.name}",
            xaxis_title="Metabolite concentration",
            yaxis_title="Residuals"
        )
        residuals_plot.add_hline(y=0, line_width=1, line_dash="dash",
                                 line_color="red")
        self._residual_plot = residuals_plot
        return self._residual_plot

    @property
    def limits(self) -> namedtuple:
        """
        Lower and upper limits for the calibration curve
        :return: namedtuple containing limits
        """

        Limits = namedtuple("Limits", "lower upper")
        return Limits(self.y.min(), self.y.max())


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
                "IS_Signal": None,
                "IS_Concentration": None
            }
        )
    )
    # Case 2: No IS, calibration points
    test_case_2 = (
        "dATP",
        pd.DataFrame.from_dict(
            {
                "Cal_name": ["SMS_01", "SMS_02", "SMS_04", "SMS_08", "SMS_16",
                             "SMS_31"],
                "Cal_Signal": np.array(
                    [89045904, 83217648, 70279200, 65918260, 37007596,
                     17540828]),
                "Cal_Concentration": np.array(
                    [2.618, 1.309, 0.655, 0.327, 0.164, 0.082]),
                "IS_Signal": None,
                "IS_Concentration": None
            }
        )
    )
    # Case 3: IS present, No [IS], no cal points
    test_case_3 = (
        "ATP",
        pd.DataFrame.from_dict(
            {
                "Cal_name": ["SMS_01", "SMS_02", "SMS_04", "SMS_08", "SMS_16",
                             "SMS_31"],
                "Cal_Signal": None,
                "Cal_Concentration": None,
                "IS_Signal": np.array(
                    [483495840, 436619616, 412101248, 381654080, 345984544,
                     316497536]),
                "IS_Concentration": None
            }
        )
    )
    # Case 4: IS present, no [IS], calibration points
    test_case_4 = (
        "ADP",
        pd.DataFrame.from_dict(
            {
                "Cal_name": ["SMS_01", "SMS_02", "SMS_04", "SMS_08", "SMS_16",
                             "SMS_31"],
                "Cal_Signal": np.array(
                    [1163372160, 506132736, 234740656, 107778792, 48474812,
                     22390422]),
                "Cal_Concentration": np.array(
                    [2.618, 1.309, 0.655, 0.327, 0.164, 0.082]),
                "IS_Signal": np.array(
                    [483495840, 436619616, 412101248, 381654080, 345984544,
                     316497536]),
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
                "IS_Signal": np.array([0.703773953]),
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

#    quant = Quantifier(test_case_4[0], test_case_4[1], calib_weight="1/X²")
#   print(quant.calibrator.equation)
# quant.calibrator.degree = 1
# print(quant.calibrator.scaled_polynome)
# # print(list(quant.calibrator.scaled_polynome.convert().coef).reverse())
#  quant.calibrator.polynome_plot.show()

# ratio = 506132736 / 436619616
# print(ratio)
# print(quant.calibrator.equation(ratio))
# print(quant.quantify(506132736, 436619616))
