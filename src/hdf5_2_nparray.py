"""Module for converting HDF5 datasets to NumPy arrays.

This module provides functionality to read data from HDF5 files into
NumPy arrays, with support for:

- Single and compound datasets
- Type conversion
- Column-wise extraction for compound data
- Memory-efficient file handling

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
May 29, 2022
"""

import h5py
import numpy as np


def hdf5_2_nparray(filename, field, dtype):
    """Read HDF5 dataset into NumPy array with type conversion.

    This function reads a specified field from an HDF5 file and converts
    it to a NumPy array. It handles both simple datasets and compound
    datasets (with multiple columns).

    Parameters
    ----------
    filename : str
        Full path to the HDF5 file
    field : str
        Name of the dataset to read from the HDF5 file.
        Can be a path within the HDF5 hierarchy (e.g. 'group/dataset')
    dtype : str or numpy.dtype
        Data type for the output array. All data will be converted
        to this type.

    Returns
    -------
    numpy.ndarray
        For simple datasets:
            1D or 2D array of specified dtype
        For compound datasets:
            2D array where each column is a field from the dataset

    Notes
    -----
    For compound datasets:
    1. Each field is extracted separately
    2. All fields are converted to the specified dtype
    3. Fields are stacked as columns in the output array
    """

    with h5py.File(filename, 'r') as hdf5_f:
        ds = hdf5_f[field]

        if ds.dtype.names is None:
            return np.asarray(ds, dtype)

        col_names = ds.dtype.names
        a = [ds[col].astype(dtype) for col in col_names]

    return np.stack(a, axis=1)
