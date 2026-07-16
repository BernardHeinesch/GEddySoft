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
    - does not follow the strategy of InnFlux (with derivatives), which was found to be inadequate with low SNR

    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
    """

    if positive_peak_only:
        b = np.hamming(filter_length)/sum(np.hamming(filter_length))
        a = 1
        filtered_cov = filtfilt(b, a, covariance.astype(float), axis=0, padtype='odd', padlen=3*(max(len(b), len(a))-1))
        idx = np.argmax(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        lag_samples = window_center - inner_window_size + idx
    else:
        if filter_length != 1:
            filtered_cov = filtfilt(np.hamming(filter_length)/sum(np.hamming(filter_length)), 1, covariance.astype(float))
        else:
            # no filtering
            filtered_cov = covariance

        if np.mean(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)]) > 0:
            # covariance is positive
            idx = np.argmax(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        else:
            # covariance is negative
            idx = np.argmin(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])

        lag_samples = window_center - inner_window_size + idx

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
