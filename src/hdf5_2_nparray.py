
import h5py
import numpy as np


def hdf5_2_nparray(filename, field, dtype):
    """
    Reads a given field of an hdf5 file into a numpy array.

    Parameters
    ----------
    filename : string
        A single string with the full hdf5 file path and name.
    field : string
        A single string with the name of the field to be extracted from the hdf5 structure.
        - If "extrapolate", then points outside the data range will be
          extrapolated.
    dtype : string
        A single string with dtype of the homogeneous output numpy array.

    comments
    --------
    Written by B. Heinesch on Sun May 29 12:00 2022
    University of Liege, Gembloux Agro-Bio Tech.
    """
    with h5py.File(filename, 'r') as hdf5_f:
        ds = hdf5_f[field]

        if ds.dtype.names is None:
            return np.asarray(ds, dtype)

        col_names = ds.dtype.names
        a = [ds[col].astype(dtype) for col in col_names]

    return np.stack(a, axis=1)
