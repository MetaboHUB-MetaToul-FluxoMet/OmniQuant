"""
Module that defines the Dataset class for input/output operations.
"""
import pathlib
from typing import Union

import pandas as pd


def _check_base_columns(data: pd.DataFrame):
    """
    This method is used to check if the data contains the required columns.

    :param data: The data to check.
    """
    essential_cols = ["Sample", "Sample Type", "Molecule", "Signal"]
    wrong_cols = [col for col in essential_cols if col not in data.columns]
    if wrong_cols:
        raise ValueError(
            f"These columns are missing: {wrong_cols}. Available "
            f"columns: {data.columns}"
        )


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

    def get_molecule(self, molecule):
        return self.data[self.data.Molecule == molecule]

    @property
    def data(self):
        return self._data

    @property
    def molecules(self):
        return self.data.Molecule.unique()

    @data.setter
    def data(self, value):
        if not isinstance(value, pd.DataFrame):
            raise ValueError("Data must be a pandas DataFrame")
        _check_base_columns(value)
        self._data = value

    @classmethod
    def from_file(cls, file: str):
        """
        This method is used to read data from a file and return it as a
        pandas DataFrame.

        :param file: The path to the file to read.
        :return: A DataSet object containing the data.
        """
        data = pd.read_csv(file, sep="\t")
        return cls(data) if data is not None else None

    def to_file(self, path: Union[str, pathlib.Path], title: str):
        """
        This method is used to write the data to a file.

        :param path: The path to the file to write.
        :param title: The title of the file.
        """
        file = pathlib.Path(path) / f"{title}.tsv"
        self.data.to_csv(file, sep="\t", index=False)

    def mutate(self, col):
        """
        This method is used to mutate the data by adding a new column.

        :param col: The column to add.
        """
        if isinstance(col, pd.Series):
            self.data = pd.concat([self.data, col], axis=1)
        else:
            raise ValueError("Column must be a pandas Series")
