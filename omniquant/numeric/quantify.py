import pandas as pd
import numpy as np

from .calibrate import Calibrator


class Quantifier:
    CASES = [1, 2, 3, 4, 5]

    def __init__(
            self,
            metabolite_name: str,
            cal_data: pd.DataFrame,
            calib_weight: str = None,
            **kwargs
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
        self.calib_weight = calib_weight
        self.response_factor = 1
        self.case = self.select_case()

        # Store any additional keyword arguments to pass to calibrator
        self._kwargs = kwargs

    def __repr__(self):

        return (f"Metabolite Name: {self.metabolite_name}\n"
                f"Calibration Data: {self.cal_data}\n"
                f"Is an internal standard present: {'Yes' if self.is_int_std else 'No'}\n"
                f"Is the internal standard concentration known: {'Yes' if self.is_int_std_conc_known else 'No'}\n"
                f"Are there multiple calibration points: {'Yes' if self.is_cal_points else 'No'}")

    @classmethod
    def from_calibration(cls, file: str, **kwargs):
        """
        This method is used to get the calibration data from the user. It
        reads the calibration data from a text file and returns it as a
        pandas DataFrame.

        :return: A pandas DataFrame containing the calibration data.
        """
        data = pd.read_csv(file, sep="\t")
        return Quantifier(cal_data=data, **kwargs)

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
        This method is used to determine the use case based on the calibration data provided.
        It checks for the presence of calibration points, internal standard, and known internal standard concentration.
        Based on these checks, it assigns a case number from 1 to 5.

        Case 1: No calibration points, no internal standard, no known internal standard concentration.
        Case 2: Calibration points present, no internal standard, no known internal standard concentration.
        Case 3: No calibration points, internal standard present, no known internal standard concentration.
        Case 4: Calibration points and internal standard present, no known internal standard concentration.
        Case 5: No calibration points, internal standard present, known internal standard concentration.

        :return: An integer representing the use case.
        """

        case = 1
        if not self.cal_data.empty:
            # Check if there are any calibration points and if there are more than one
            self.is_cal_points = self.cal_data["Signal"].any() and len(
                self.cal_data["Signal"] > 1)
            # Check if an internal standard is present
            self.is_int_std = self.cal_data["IS_signal"].any()
            # Check if the concentration of the internal standard is known
            self.is_int_std_conc_known = self.cal_data[
                "IS_Concentration"].any()

            # Assign the case number based on the checks
            if self.is_cal_points:
                case = 2
            if self.is_int_std:
                case = 3
                if self.is_cal_points:
                    case = 4
            if self.is_int_std_conc_known:
                case = 5
        return case

    def _set_quantifier_func(self):
        """
        This method sets the quantification function based on the selected
        use case. It uses a dictionary to map the case numbers to the
        corresponding quantification methods. The quantification function is
        then set by looking up the method corresponding to `self.case` in the
        dictionary.
        """
        func_map = {
            1: self._quantify_no_int_std_no_curve,
            2: self._quantify_no_int_std_with_curve,
            3: self._quantify_int_std_no_conc_no_curve,
            4: self._quantify_int_std_no_conc_with_curve,
            5: self._quantify_int_std_with_conc
        }
        self.quantify = func_map.get(self.case)

    @property
    def calibrator(self):
        """
        Initialize or get the calibrator. Only used in cases 2 & 4.
        :return: Initialized calibrator object
        """

        if self._calibrator:
            return self._calibrator
        if self.case not in [2, 4]:
            raise CaseError(
                f"Use case is initialized as {self.case}. Only cases 2 & 4 get a calibrator")

        try:
            self._calibrator = Calibrator(
                name=self.metabolite_name,
                x=self.cal_data["Cal_Concentration"].to_numpy(),
                y=self.cal_data["Signal"].to_numpy() if not self.is_int_std
                else np.divide(self.cal_data["Signal"].to_numpy(),
                               self.cal_data["IS_signal"]),
                case=self.case,
                weight=self.calib_weight,
                **self._kwargs
            )
        except Exception as e:
            raise CalibrationError(
                f"There was an error while initializing the calibrator: {e}")

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
    def _quantify_no_int_std_no_curve(signal, reference=1):
        """
        Relative quantification of our signal.
        :param signal: Metabolite sample signal/area
        :return: Metabolite sample signal/area
        """
        if not isinstance(signal, float) and not isinstance(signal, int):
            raise ValueError(
                f"Signal should be a number. Detected type: {type(signal)}:")
        return signal / reference

    def _quantify_no_int_std_with_curve(self, signal):
        """
        Quantify signal using a calibration curve and no internal standard
        signal.
        :param signal: sample signal
        :return: tuple(value, flag) where value is the quantified value and
        flag is a string indicating if the value is above the ULOQ or below
        the LLOQ. If the value is within the calibration range, flag is None.
        """
        return self._quantify_int_std_no_conc_with_curve(signal, 1)

    def _quantify_int_std_no_conc_no_curve(self, signal):
        """
        Quantify signal using an internal standard signal and no calibration curve.
        :param signal:
        :return:
        """
        if len(self.calibrator.Cal_data["IS_Signal"]) > 1:
            raise QuantificationError(
                "Expected only one value for the  internal signal. Number "
                "of detected  values:"
                f" {len(self.calibrator.Cal_data['IS_Signal'])}")
        return signal / self.calibrator.Cal_data["IS_Signal"]

    def _quantify_int_std_no_conc_with_curve(self, signal, std_signal):
        """
        Quantify signal using a calibration curve and internal standard signal.
        :param signal: sample signal
        :param std_signal: sample internal standard signal
        :return: tuple(value, flag) where value is the quantified value and
        flag is a string indicating if the value is above the ULOQ or below
        the LLOQ. If the value is within the calibration range, flag is None.
        """
        to_quant = signal / std_signal
        if self.calibrator.degree == 2:
            res = (self.calibrator.equation - to_quant).roots[1]
        else:
            res = (to_quant - self.calibrator.equation[0]) / \
                  self.calibrator.equation[1]

        if to_quant > self.calibrator.limits.upper:
            return res, ">ULOQ"
        if to_quant < self.calibrator.limits.lower:
            return res, "<LLOQ"
        return res, None

    @apply_response_factor
    def _quantify_int_std_with_conc(self, signal):
        """
        This method is used to quantify the signal using the known
        concentration of the internal standard. It multiplies the signal
        with the known concentration of the internal standard and divides it
        by the internal standard signal to get the quantified value.

        :param signal: The signal of the sample to be quantified.
        :return: The quantified value of the signal.
        """
        return (signal / self.cal_data["IS_Signal"]) * self.cal_data[
            "IS_Concentration"]


class CaseError(Exception):
    pass


class CalibrationError(Exception):
    pass


class QuantificationError(Exception):
    pass
