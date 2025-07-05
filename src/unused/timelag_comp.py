# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 17:15:45 2022

@author: Ariane FAURES

Functions to compute the covariance and then the timelag of HF eddy covariance data
TO DO:
    - Add the possibility to to linear detrending
    - add the hlag computation

"""
import numpy as np
import pandas as pd

# %% List of functions

"""
Fonction that computes the fluctuations around the mean of a time series using 
a given method.
Inputs:
    - data: [pd dataframe] with the data
    - detrend: [string] method (for the moment only BA)
    - c, h, w: [string] column names of the CO2, H2O and wind variables in the 
            dataframe
Outputs:
    - x_prime: [pd dataframe] with the c, h or w fluctuation
"""
def fluctuations (data,detrend,c,h,w):
    if detrend == 'BA':
        c_mean = np.mean(data[c])
        h_mean = np.mean(data[h])
        w_mean = np.mean(data[w])
        c_prime = data[c] - c_mean
        h_prime = data[h] - h_mean
        w_prime = data[w] - w_mean
    
    return(c_prime, h_prime, w_prime)


"""
Fonction that computes the covariance wc and wh as a function of a lag.
Note: for the moment wc only. Add wh
Inputs:
    - x_prime: [pd dataframe] with the c, h or w fluctuation
    - win_cen: [int] centre of the lag window
    - win_min, win_max : [int] min and max of the total lag window
    - win_halfsize: [int] half size of the window in which the lag will be 
            searched for. So this window will be centre +- halfsize
Outputs:
    - cov_lag: [pd dataframe] covariance as a function of the lag in samples
    - cov_max: [list] extremum of the covariance function in the given window
            with the corresponding samples lag
"""
def timelag_comp(c_prime,h_prime,w_prime,win_cen,win_min, win_max, win_halfsize):
    cov_lag = pd.DataFrame({'Samples': [0], 'cov':[0]})
    k=0
    for i in [*range(win_min,win_max+1)]: 
        if i == 0:
            cov = np.mean(c_prime.to_numpy()*w_prime.to_numpy())
            cov_lag.loc[k] = [i,cov]
            k = k+ 1
        elif i > 0:
            cov = np.mean(c_prime.to_numpy()[i:]*w_prime.to_numpy()[0:-i])
            cov_lag.loc[k] = [i,cov]
            k = k+ 1
        elif i < 0:
            cov = np.mean(w_prime.to_numpy()[-i:]*c_prime.to_numpy()[0:i])
            cov_lag.loc[k] = [i,cov]
            k = k+ 1
    
    cov_lag = cov_lag.set_index('Samples') # Set the lag samples as the dataframe index
    cov_lag_abs = abs(cov_lag['cov'])
    cov_max = cov_lag_abs.loc[win_cen-win_halfsize:win_cen+win_halfsize].max()
    idx_max = cov_lag_abs.loc[win_cen-win_halfsize:win_cen+win_halfsize].idxmax()
    cov_max = [idx_max,cov_max]
    #in_win = np.arange(win_cen-win_halfsize,win_cen+win_halfsize,1)
    
    return(cov_lag,cov_max)
