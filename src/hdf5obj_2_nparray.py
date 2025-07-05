import numpy as np


def hdf5obj_2_nparray(hdf5obj, dtype):
    """
    Converts a given hdf5 dataset object into a numpy array.

    Parameters
    ----------
    hdf5obj : Dataset object
        Dataset object to be converted to numpy array
    dtype : string
        A single string with dtype of the homogeneous output numpy array.

    comments
    --------
    Written by B. Heinesch on Sun May 29 12:00 2022
    University of Liege, Gembloux Agro-Bio Tech.
    """

    if hdf5obj.dtype.names is None:
        return np.asarray(hdf5obj, dtype)

    col_names = hdf5obj.dtype.names
    a = [hdf5obj[col].astype(dtype) for col in col_names]

    return np.stack(a, axis=1)
