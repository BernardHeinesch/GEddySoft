# %%
import numpy as np
import time
import csv

def xcov(x,y,lag):
    """
    Perform Cross-Correlation on x and y
    x    : 1st signal
    y    : 2nd signal
    lag  : lag interval [start,end] in samples

    returns
    corr : lagged covariance
    
    comments:
        - when translating from MATLAB (InnFlux), no xcov function was found in python packages 
        - MATLAB computes the mean of the initial array while python (and excel) uses the mean 
          of the truncated arrays. So gives similar results for a zero lag but different results 
          for a non-zero lag !
        
    """

    lag = np.linspace(lag[0], lag[1], num=lag[1]-lag[0]+1, dtype='int16')

    # create output array
    crosscov = np.full((len(lag)), np.nan)

    for i in range(0, len(lag), 1):

        yshifted = np.roll(y, lag[i])
        if lag[i] > 0:
            crosscov[i] = ((x[lag[i]:] * yshifted[lag[i]:]).sum() - x[lag[i]:].sum() * yshifted[lag[i]:].sum() / x[lag[i]:].size) / (x[lag[i]:].size - 1)
        else:
            crosscov[i] = ((x[0:len(x)-1-lag[i]] * yshifted[0:len(x)-1-lag[i]]).sum() - x[0:len(x)-1-lag[i]].sum() * yshifted[0:len(x)-1-lag[i]].sum() / x[0:len(x)-1-lag[i]].size) / (x[0:len(x)-1-lag[i]].size - 1)

    # reduce to float if unique lag sent
    if len(crosscov) == 1:
        crosscov = crosscov[0]

    return crosscov

import numpy as np
from scipy.signal import filtfilt
import matplotlib.pyplot as plt


def find_covariance_peak(covariance, outer_window_size, inner_window_size, window_center, filter_length, positive_peak_only, plot=False, filename=''):
    """
    For two shifted (anti)correlated signals, s1 and s2,
    this fuction estimates the lag of s2 with respect to s1 by finding the
    maximum/minimum of the (smoothed) crosscovariance xcov(s1, s2)

    parameters
    ----------
    covariance: np, covariance function
    outer_window_size: samples, size (+/-) of the outer window (for lag plotting only)
    inner_window_size: samples, size (+/-) of the inner window (for lag search)
    window_center: samples, lag of the outer window center
    filter_length: samples smoothing length of filter applied prior to finding covariance peak
    positive_peak_only: boolean, True if only the positive peak in the cov function has to be searched for

    returns
    -------
    lag_samples: lag in samples

    Comments
    --------
    - if filter_lenght is too high, an error will occur in the filtfilt function
    - code in comments is the original one from INNFLUX but was shown inappropriate for low signal-to-noise ratio

    # """

    if positive_peak_only:
        b = np.hamming(filter_length)/sum(np.hamming(filter_length))
        a = 1
        filtered_cov = filtfilt(b, a, covariance.astype(float), axis=0, padtype='odd', padlen=3*(max(len(b), len(a))-1))
        idx = np.argmax(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        lag_samples = window_center - inner_window_size + idx + 1
        # BH : check if +1 is needed
    else:
        # # find maximum or minimum of filtered derivative
        # b = np.ones(100, np.int8)/100
        # a = 1
        # temp = np.diff(covariance, n=1, axis=0).astype(float)
        # padlen = 3*(max(len(b), 1)-1)
        # if padlen >= len(temp): padlen = len(temp)-1
        # dfdt = filtfilt(b, a, temp, axis=0, padtype='odd', padlen=padlen)
        # dfdt_min_idx = np.argmin(dfdt[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        # dfdt_max_idx = np.argmax(dfdt[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])

        if filter_length != 1:
            filtered_cov = filtfilt(np.hamming(filter_length)/sum(np.hamming(filter_length)), 1, covariance.astype(float))
        else:
            # no filtering
            filtered_cov = covariance

        # # find the index of the extremum (max or min depending on the concavity of the covariance function)
        # if dfdt_max_idx < dfdt_min_idx:
        #     # maximum of derivative is to the left -> covariance is positive
        #     idx = np.argmax(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        # else:
        #     # maximum of derivative is to the right -> covariance is negative
        #     idx = np.argmin(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])


        if np.mean(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)]) > 0:
            # covariance is positive
            idx = np.argmax(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        else:
            # covariance is negative
            idx = np.argmin(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])


        lag_samples = window_center - inner_window_size + idx + 1

    # plot results if requested
    if plot:

        fig, ax = plt.subplots(2, 1, figsize=(8, 6))

        # plot whole cov function and found maximum
        xrange = np.linspace(-outer_window_size, outer_window_size, num=outer_window_size*2+1) + window_center
        ax[0].plot(xrange, covariance)
        ax[0].plot(xrange, filtered_cov)
        ax[0].set_ylabel('Cov (uncal units)')
        ax[0].set_xlabel('Samples')

        ax[0].set_title('Lag determination for file :' + filename)

        ax[0].axvline(x=window_center, color='black', label='axvline - full height')
        ax[0].axvline(x=-inner_window_size+window_center, color='black', linestyle='dashed')
        ax[0].axvline(x=inner_window_size + window_center, color='black', linestyle='dashed')
        ax[0].plot(lag_samples, filtered_cov[outer_window_size + lag_samples - window_center], marker="o", markersize=8, markeredgecolor="red", markerfacecolor="red")

        # zoom and search range
        xrange = np.linspace(-inner_window_size, inner_window_size, num=(2*inner_window_size)) + window_center
        ax[1].plot(xrange, covariance[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        ax[1].plot(xrange, filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        ax[1].set_ylabel('Cov (uncal units)')
        ax[1].set_xlabel('Samples')

        ax[1].axvline(x=window_center, color='black')
        ax[1].axvline(x=-inner_window_size + window_center, color='black', linestyle='dashed')
        ax[1].axvline(x=inner_window_size + window_center, color='black', linestyle='dashed')
        ax[1].plot(lag_samples, filtered_cov[outer_window_size + lag_samples - window_center], marker="o", markersize=8, markeredgecolor="red", markerfacecolor="red")

        fig.subplots_adjust(hspace=0.3)

        plt.show(block=False)
        # plt.pause(1)  # pauses one second before disappearing
        # plt.close('all')  # close all previous plots

    return int(lag_samples)

"""
Read the data and metadata files extracted with the unzip function  
Inputs:
    - unzipped data file name (from unzip_GHG function)
    - unzipped metadata file name (")
Outputs:
    - pandas dataframe with high frequency data
    - pandas dataframe with header of data file
    - pandas dataframe with metadata information
"""
import pandas as pd

def read_GHG (data_name,metadata_name): 
        # Read the header (first 7 lines of the file)
        # And afterwards read the many body (with variable names as columns)
        with open(data_name, mode='r') as file:
            file_header = pd.read_table(file, nrows = 6, header = None)
        with open(data_name, mode='r') as file:    
            file_data = pd.read_table(file, header = 7)
            
        file.close()
            
        return([file_header,file_data])


file_header,file_data = read_GHG(r'D:\z-DATA\Database\Lonzee\2_Eddy covariance\3_LICOR-7200\0_GHG\2023-01-01T050000_OTL.data',r'D:\z-DATA\Database\Lonzee\2_Eddy covariance\3_LICOR-7200\0_GHG\2023-01-01T050000_OTL.metadata')

# the extremum of the lag function (smoothed according to LAG_COVPEAK_FILTER_LENGTH parameter)
# will be search for, in the lag window [LAG_EXPECTED - LAG_INNER_WINDOW_SIZE, LAG_EXPECTED + LAG_INNER_WINDOW_SIZE].
# if an extremum of the covariance is not found inside the window (but at one of its extremes), the time lag is set to LAG_EXPECTED.

# expected physical lag (in samples of sampling_rate_final)
LAG_TRACER_WINDOW_CENTER = -4
# time lag plotting window around the expected lag (in samples of sampling_rate_final)
# ex: if = 100, the plotting window will be [LAG_EXPECTED-100, LAG_EXPECTED+100]
LAG_OUTER_WINDOW_SIZE = 100
# time lag search window around the expected lag (in samples of sampling_rate_final)
# ex: if = 20, the search window will be [LAG_EXPECTED-20, LAG_EXPECTED+20]
LAG_INNER_WINDOW_SIZE = 20
# smoothing length of filter applied prior to finding covariance peak, in samples of sampling_rate_final
LAG_COVPEAK_FILTER_LENGTH = 1
# plotting option (1 for plotting)
PLOT_FIND_COVARIANCE_PEAK = 1

# calculate cross-covariance over a given lag-window
cov_wc = xcov(file_data['W (m/s)'],file_data['CO2 dry(umol/mol)'], [LAG_TRACER_WINDOW_CENTER - LAG_OUTER_WINDOW_SIZE,LAG_TRACER_WINDOW_CENTER + LAG_OUTER_WINDOW_SIZE])

# find time lag
lag_samples = find_covariance_peak(cov_wc, LAG_OUTER_WINDOW_SIZE, LAG_INNER_WINDOW_SIZE, LAG_TRACER_WINDOW_CENTER, LAG_COVPEAK_FILTER_LENGTH, 0, PLOT_FIND_COVARIANCE_PEAK)  # lag in samples
# if lag was found at the extremum of the search window, replace by the default
if (lag_samples == LAG_TRACER_WINDOW_CENTER - LAG_INNER_WINDOW_SIZE + 1 or 
    lag_samples == LAG_TRACER_WINDOW_CENTER + LAG_INNER_WINDOW_SIZE):
    lag_samples = LAG_TRACER_WINDOW_CENTER

# shift concnetrations to correct for the lag
c_prime_shifted = np.roll(file_data['CO2 dry(umol/mol)'], lag_samples)


