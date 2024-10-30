"""
Module that defines the Dataset class for input/output operations.
"""
import pathlib
from typing import Union

import pandas as pd


class DataSet:
    """
    Class that defines the DataSet object for input/output operations.
    """

    def __init__(self, data: pd.DataFrame):
        self.data = data

    def __repr__(self):
        return self.data.__repr__()

    def __call__(self):
        """
        On call return the pandas dataframe
        :return: pd.DataFrame containing data
        """
        return self.data

    def get_molecule(self, analyte: str):
        return self.data[self.data.Analyte == analyte]

    @property
    def data(self):
        return self._data

    @property
    def analytes(self):
        return self.data.Analyte.unique()

    @data.setter
    def data(self, value):
        if not isinstance(value, pd.DataFrame):
            raise ValueError("Data must be a pandas DataFrame")
        self._data = value

    @classmethod
    def from_file(cls, path: str) -> "DataSet":
        """
        This method is used to read data from a file and return it as a
        pandas DataFrame.

        :param path: The path to the file to read.
        :return: A DataSet object containing the data.
        """
        data = pd.read_csv(path, sep="\t")
        return cls(data) if data is not None else None

    def to_file(self, path: Union[str, pathlib.Path], title: str):
        """
        This method is used to write the data to a file.

        :param path: The path to the file to write.
        :param title: The title of the file.
        """
        file = pathlib.Path(path) / f"{title}.tsv"
        self.data.to_csv(file, sep="\t", index=False)

    # def mutate(self, col):
    #     """
    #     This method is used to mutate the data by adding a new column.
    #
    #     :param col: The column to add.
    #     """
    #     if isinstance(col, pd.Series):
    #         self.data = pd.concat([self.data, col], axis=1)
    #     else:
    #         raise ValueError("Column must be a pandas Series")


class SampleData(DataSet):

    def __init__(self, data: pd.DataFrame):
        super().__init__(data)
        self._check_sample_data()

    def _check_sample_data(self):
        """
        This method is used to check if the data contains the required columns.
        """
        essential_cols = ["Sample", "Sample Type", "Analyte", "Signal"]
        wrong_cols = [col for col in essential_cols if col not in
                      self.data.columns]
        if wrong_cols:
            raise ValueError(
                f"These columns are missing from sample data: {wrong_cols}. "
                f"Available columns: {self.data.columns}"
            )

    @property
    def cal_points(self):
        return self.data[self.data["Sample Type"].str.contains("Cal")]

    @cal_points.setter
    def cal_points(self, value: pd.DataFrame):
        self.data[self.data["Sample Type"].str.contains("Cal")] = pd.concat(
            objs=[self.data[self.data["Sample Type"].str.contains("Cal")],
                  value],
            axis=0
        )

    @property
    def sample_points(self) -> pd.DataFrame:
        return self.data[self.data["Sample Type"] == "Unknown"]

    @sample_points.setter
    def sample_points(self, value: pd.DataFrame):
        self.data[self.data["Sample Type"] == "Unknown"] = pd.concat(
            objs=[self.data[self.data["Sample Type"] == "Unknown"], value],
            axis=0
        )


class CalibrationData(DataSet):

    def __init__(self, data: pd.DataFrame):
        super().__init__(data)
        self.check_calibration_data()
        self.response_factors = self.isolate_response_factors()

    def check_calibration_data(self):
        """
        This method is used to check if the calibration data is valid.
        """
        if "Analyte" not in self.data.columns:
            raise ValueError("The calibration data must contain the "
                             "'Analyte' column.")

    def isolate_response_factors(self) -> Union[dict, None]:
        """Isolate the response factors from the calibration data."""

        if "Response Factor" in self.data.columns:
            return self.data["Response Factor"].to_dict()
        return None
