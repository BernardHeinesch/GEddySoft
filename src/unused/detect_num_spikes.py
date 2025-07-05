import pandas as pd
import numpy as np

def detect_num_spikes(x, threshold_stds, window_width):
    """
    Spikes detection
    Consider improving this routine by using the VM97 algorithm, described and coded in R in Vitale et al. 2020, Biogeosciences, "Robust data cleaning procedure
    for eddy covariance flux measurements"", routine "despiking"

    parameters
    ----------
    x :  raw high frequency eddy covariance time series.
    threshold_stds : number of times the std dev to define a spike
    window_width : size of the rolling windows for medion and std computation (in samples)

    returns
    -------
    num_spikes: the number of detected spikes
    
    Comments
    --------

    Written by B. Heinesch, 10 April, 2023.
    University of Liege, Gembloux Agro-Bio Tech.

    """
    xdf = pd.DataFrame(x)    
    xmedian = xdf.rolling(window=window_width).median()
    xstd = xdf.rolling(window=window_width).std()

    n=abs(xdf - xmedian) > threshold_stds*xstd

    num_spikes = n.values.sum()
    
    return num_spikes
