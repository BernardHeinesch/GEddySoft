# -*- coding: utf-8 -*-
"""
Created on Mon May 22 17:29:24 2023

@author: u230638
"""

import math
import xarray as xr
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
from datetime import datetime
from datetime import timedelta
import configparser
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
import re
import sys
from collections import defaultdict
import seaborn as sns
from matplotlib.patches import Patch
import h5py

from TOFEC_read_results import read_results
from TOFEC_plot_xarray import plot_ds_time, plot_ds_corr, plot_ds_cov
from TOFEC_plot_xarray import multiplot_ds_diel_period, plot_ds_diel_period
from TOFEC_plot_xarray import plot_ds_clock, plot_ratio_mz
from TOFEC_func import mz_attr, verif_concat, mean_xh, multiplot
from TOFEC_basic_QC import visual_timelag, visual_stationnarity, keep_stat_fluxes
from TOFEC_basic_QC import visual_raw_fluxes
from TOFEC_signif_mz import strat_Langford_LOD, strat_mean_LOD, strat_short_exchange, save_ds_sign
from TOFEC_correlogram import flux_correlogram
from TOFEC_calib import apply_calib_Crist, apply_calib_kinetic, atom_concentration
from TOFEC_plumes_remove import remove_plumes
from TOFEC_budget import budget_bar_plot, budget_stacked_plot, pie_plot, relative_contribution, bvoc_fingerprint
from TOFEC_chemistry import kendrick_plot
from TOFEC_list_files import list_files_datestring
from TOFEC_concat_hdf5 import concat_hdf5_new

sns.set(color_codes=True)

# =============================================================================
# %% Read code config
# =============================================================================
config = configparser.ConfigParser()
config.read('BE-Lon_N2O_EC_config.ini')

# =============================================================================
# %% Verify concatenation
# If GEddySoft daily output files have not already been concatenated,
# run TOFEC_concat_hdf5.py
# =============================================================================
verif_concat(config)

# =============================================================================
# %% Read VOC, O3 and ICOS fluxes files
# Import all into ds + ds_cov
# =============================================================================
print('-----------------------------------------------------------------------------')
print('Read VOC, ICOS, O3 files')
print('-----------------------------------------------------------------------------')

ds = read_results(config)[0]

# # Add LOD computation
# ds['LOD'] = abs(int(config['ANALYSIS_OPTIONS']['LOD_mult_factor'])*ds['random_error_FS'])
# ds['LOD'].attrs['name'] = 'flux limits of detection'
# ds['LOD'].attrs['units'] = 'ug m-2 s-1'

# # Add concentration variance
# ds['conc_var'] = ds['conc_std']**2
# ds['conc_var'].attrs['units'] = 'ppbv² or uncal'
# ds['conc_var'].attrs['name'] = 'concentration variance'

# # Remove HH affected by plumes
# if config['ANALYSIS_OPTIONS']['plumes_remove'] == 'True':
#     ds = remove_plumes(ds, config, option_removal='var_crit')


if config['USER_OPTIONS']['concat_cov'] == 'True':
    ds_cov = read_results(config)[1]

# mz_attr_list = ['mz_glov', 'mz_match', 'mz_chemform', 'mz_name']
# mz_names = mz_attr(ds, 'mz_chemform', 'all')

# # Compute deposition velocity
# ds['flux_molar'] = ds['flux']*1000./(ds.mz-1)
# ds['flux_molar'].attrs['name'] = 'Flux in molar or uncal basis'
# ds['flux_molar'].attrs['units'] = 'nmol m-2 s-1 or ncps mol m-2 s-1'

# # Use equation from EddyPro /!\ use the sonic temperature as Bernard
# ds['ambient_molar_vol'] = 8.314*ds['T']/101325.
# ds['ambient_molar_vol'].attrs['name'] = 'Average ambient air molar volume'
# ds['ambient_molar_vol'].attrs['units'] = 'm3 mol-1'

# ds['flux_cinematic'] = ds['flux_molar']*ds['ambient_molar_vol']
# ds['flux_cinematic'].attrs['name'] = 'Cinematic flux'
# ds['flux_cinematic'].attrs['units'] = 'ppbv m s-1 or ncps m s-1'

# # When computing deposition velocitiy, avoid dividing by too low concentrations
# ds['dep_velocity'] = ds['flux_cinematic']/ds['conc_mean'].where(
#     ds['conc_mean'] >= ds['conc_mean'].quantile(0.1, dim='time'))
# # in cm s-1
# ds['dep_velocity'] = ds['dep_velocity']*-100.
# ds['dep_velocity'].attrs['name'] = 'Deposition velocity'
# ds['dep_velocity'].attrs['units'] = 'cm s-1'

# =============================================================================
# %% Read Fluxnet data
# =============================================================================
os.chdir(config['FILES_PATHS']['L2_O3_folder'])
# Read L2 file
df_fluxnet = pd.read_csv('ICOSETC_BE-Vie_FLUXNET_INTERIM_HH_L2' + '.csv')
df_fluxnet.replace(to_replace=-9999, value=np.nan, inplace=True)
df_fluxnet.index = pd.to_datetime(df_fluxnet['TIMESTAMP_END'],
                             format='%Y%m%d%H%M').rename('time')
# Select time values in df_fluxnet that are similar to ds
df_fluxnet = df_fluxnet.merge(ds.time.to_dataframe(), left_index=True, right_index=True)

# List of variables from df_L2 to add into ds
list_var = pd.read_csv('ICOSETC_BE-Vie_FLUXNET_INTERIM_HH_L2_list_var'+'.csv')

for var, new_var, name, unit in zip(list_var.Var, list_var.New_var, list_var.Name, list_var.Unit):
    ds[new_var] = df_fluxnet[var].to_xarray()
    ds[new_var].attrs['name'] = name
    ds[new_var].attrs['units'] = unit
# =============================================================================
# %% Remove HH in the rear of the sonic anemometer
# =============================================================================
if config['ANALYSIS_OPTIONS']['sonic_filtering'] == 'True':
    print('-----------------------------------------------------------------------------')
    print('Remove HH in the rear of the sonic anemometer')
    print('-----------------------------------------------------------------------------')
    rear_anemometer = float(config['ANALYSIS_OPTIONS']['north_offset']) - 180
    range_filtering = float(config['ANALYSIS_OPTIONS']['range_wind_filtering'])
    mask_filtering = (ds.WD_1_1_1 < (rear_anemometer - range_filtering))\
        | (ds.WD_1_1_1 > (rear_anemometer + range_filtering))
    # Filter ds
    ds = ds.where(mask_filtering)

# =============================================================================
# %% For BH: data visualization
# =============================================================================
# Select variable, mz, period
var_plot = 'flux'
var_mz = True # mz dependent variable (True) or undependent (False)
mz_plot = 69.07
tstart = '2023-05-01' # e.g. '2022-05-01'
tend = '2023-06-01' # e.g. '2022-06-01'

# Select subset of dataset ds
ds_sel = ds.sel(time=slice(tstart, tend))

# Plot
plt.figure()
if var_mz is True:
    plt.scatter(ds_sel.time, ds_sel.sel(mz=mz_plot, method='nearest')[var_plot])
else:
    plt.scatter(ds_sel.time, ds_sel[var_plot])
plt.ylabel(ds_sel[var_plot].attrs['name'] + '\n' + ds_sel[var_plot].attrs['units'])

# =============================================================================
# %% Remove plume events data
# =============================================================================
wind_inf = 180
wind_sup = 300

# Remove data flagged by the criterion
if config['ANALYSIS_OPTIONS']['plumes_remove'] == 'True':
    # Compute 99th percentile for "non-polluted" wind sectors
    ds['q99'] = ds['conc_var'].where((ds.WD_1_1_1 < wind_inf) |
                                     (ds.WD_1_1_1 > wind_sup)).groupby('mz').quantile(q=0.99, dim='time')
    # Create a flag when conc_var is higher than q99 in the "polluted" wind sector
    ds['flag_q99'] = ds['conc_var']*0.
    ds['flag_q99'] = xr.where((ds['conc_var'] > ds['q99'])\
                              & (ds.WD_1_1_1 > wind_inf)\
                                  & (ds.WD_1_1_1 < wind_sup),
                                  1, 0)
    print('-----------------------------------------------------------------------------')
    print('Remove plume events')
    print('-----------------------------------------------------------------------------')
    for var in ['conc_mean', 'flux', 'dep_velocity']:
        ds[var] = ds[var].where(ds['flag_q99'] != 1, np.nan)

# =============================================================================
# %% Remove outliers
# =============================================================================
if config['ANALYSIS_OPTIONS']['remove_outliers'] == 'True':
    print('-----------------------------------------------------------------------------')
    print('Remove outliers')
    print('-----------------------------------------------------------------------------')
    
    ds['flux'] = ds['flux'].where(
        (ds['flux'] >= ds['flux'].quantile(0.001, dim='time')) &
        (ds['flux'] <= ds['flux'].quantile(0.999, dim='time')),
        np.nan)

# =============================================================================
# %% Fluxes calibration
# =============================================================================
if config['ANALYSIS_OPTIONS']['apply_calib'] == 'True':
    print('-----------------------------------------------------------------------------')
    print('Fluxes calibration')
    print('-----------------------------------------------------------------------------')
    ds_cal = apply_calib_Crist(ds, config)
    ds_cal_kinetic = apply_calib_kinetic(ds)

# =============================================================================
# %% Basic QC
# =============================================================================
if config['ANALYSIS_OPTIONS']['basic_QC'] == 'True':
    print('-----------------------------------------------------------------------------')
    print('Quality check')
    print('-----------------------------------------------------------------------------')
    # Create a new folder to store results
    if not os.path.exists(config['FILES_PATHS']['output_results_folder']+'\\QC'):
        os.makedirs(config['FILES_PATHS']['output_results_folder']+'\\QC')
    # =============================================================================
    # Timelag
    # =============================================================================
    ds['lagtime_phys'] = ds['lagtime']-ds['lagtime_instr']
    visual_timelag(ds, config)
    # =============================================================================
    # Stationnarity
    # =============================================================================
    # Plot stationnarity
    visual_stationnarity(ds, config, 'SST_FW96')
    if config['ANALYSIS_OPTIONS']['keep_stat'] == 'True':
        ds = keep_stat_fluxes(ds, config, 'SST_FW96')
    # =============================================================================
    # Raw fluxes
    # =============================================================================
    visual_raw_fluxes(ds, config)

# =============================================================================
# %% Covariance functions multiplot
# =============================================================================
if config['ANALYSIS_OPTIONS']['cov_func'] == 'True':
    print('-----------------------------------------------------------------------------')
    print('Covariance functions')
    print('-----------------------------------------------------------------------------')
    output_folder = config['FILES_PATHS']['output_results_folder']+'\\Cov_functions'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Chose start and end dates between which cov_func must be plotted
    start_date = datetime.strptime(eval(config['ANALYSIS_OPTIONS']['date_cov_func'])[0],
                                   "%Y_%m_%d__%H_%M_%S")
    end_date = datetime.strptime(eval(config['ANALYSIS_OPTIONS']['date_cov_func'])[1],
                                 "%Y_%m_%d__%H_%M_%S")
    # Create HH dates array
    date_array = []
    current_date = start_date
    while current_date <= end_date:
        date_array.append(current_date)
        current_date += timedelta(minutes=30)

    # Loop and plot covariances functions
    for date in date_array:
        print(date)
        if np.all(np.isnan(ds_cov.sel(time=date, method='nearest').cov_func)).values:
            print('No covariance function')
        else:
            print('Covariance function plotted')
            # Plot of covariance functions
            plot_ds_cov(ds_cov=ds_cov,
                    mz_list='all',
                    HH_cov=date,
                    concat_option=config['USER_OPTIONS']['concat_option'],
                    suptitle='BVOC covariance functions',
                    mz_info='mz_chemform',
                    kind='scatter',
                    numerical_info=None,
                    figclose=True,
                    figwidth=max(20, ds_cov.mz.size*2/3),
                    figsave=True,
                    output_folder=output_folder)

# =============================================================================
# %% Important variables multiplot
# =============================================================================
if config['ANALYSIS_OPTIONS']['create_multiplot'] == 'True':
    print('-----------------------------------------------------------------------------')
    print('Fluxes multiplot')
    print('-----------------------------------------------------------------------------')
    output_folder = config['FILES_PATHS']['output_results_folder']+'\\All_mz_multiplot'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    multiplot(ds, config, output_folder)

# =============================================================================
# %% Significant fluxes
# =============================================================================
if config['ANALYSIS_OPTIONS']['apply_sign_mz'] == 'True':
    print('-----------------------------------------------------------------------------')
    print('Selection of significant m/z')
    print('-----------------------------------------------------------------------------')
    # Create a new folder to store results
    if not os.path.exists(config['FILES_PATHS']['output_results_folder']+'\\Significant'):
        os.makedirs(config['FILES_PATHS']['output_results_folder']+'\\Significant')
    if config['ANALYSIS_OPTIONS']['signif_strat'] == 'Langford_LOD':
        significant = strat_Langford_LOD(ds, config,
                                         view_flux_LODs=False,
                                         error_stat='random_error_FS',
                                         mult_factor=config['ANALYSIS_OPTIONS']['LOD_mult_factor'],
                                         save_stat_fig=False,
                                         save_stat_xlsx=True,
                                         save_sign_plot=False)
        stat_flux_S_N = significant[0]
        ds_sign = significant[1]
    if config['ANALYSIS_OPTIONS']['signif_strat'] == 'mean_LOD':
        significant = strat_mean_LOD(ds, config,
                                         view_flux_LODs=True,
                                         error_stat='random_error_FS',
                                         mult_factor=config['ANALYSIS_OPTIONS']['LOD_mult_factor'],
                                         flux_S_N_lim=0.7,
                                         save_stat_fig=True,
                                         save_stat_xlsx=True,
                                         save_sign_plot=True)
        stat_flux_S_N = significant[0]
        ds_sign = significant[1]
    elif config['ANALYSIS_OPTIONS']['signif_strat'] == 'short_exchange':
        significant = strat_short_exchange(ds, config,
                                           error_stat='random_error_FS',
                                           mult_factor=config['ANALYSIS_OPTIONS']['LOD_mult_factor'],
                                           hour_start=8,
                                           hour_end=20,
                                           flux_day_flag_ratio_lim=0.25,
                                           flux_period_size=3,
                                           save_stat_fig=True,
                                           save_stat_xlsx=True,
                                           save_sign_plot=True)
        stat_flux_S_N = significant[0]
        ds_sign = significant[1]
    
    # Save ds_sign as csv
    if config['USER_OPTIONS']['output_csv'] == 'True':
        save_ds_sign(ds_sign, config)

# =============================================================================
# %% Flux correlogram
# =============================================================================
if config['ANALYSIS_OPTIONS']['apply_sign_mz'] == 'True'\
    and config['ANALYSIS_OPTIONS']['flux_correlogram'] == 'True':
        print('-----------------------------------------------------------------------------')
        print('Flux correlogram')
        print('-----------------------------------------------------------------------------')
        # Create a new folder to store results
        if not os.path.exists(config['FILES_PATHS']['output_results_folder']+'\\Flux_correlogram'):
            os.makedirs(config['FILES_PATHS']['output_results_folder']+'\\Flux_correlogram')
        flux_correlogram(ds_sign, config)
# =============================================================================
# %% Budget bar plot
# =============================================================================
hour_start = 0
hour_end = 24

# tstart_1 = '2023-04-15'
# tend_1 = '2023-11-01'

# tstart_2 = '2023-04-15'
# tend_2 = '2023-05-15'

# tstart_3 = '2023-06-05'
# tend_3 = '2023-06-25'

# tstart_4 = '2023-07-15'
# tend_4 = '2023-08-05'

# tstart_5 = '2023-08-05'
# tend_5 = '2023-08-25'

# tstart_6 = '2023-10-15'
# tend_6 = '2023-11-01'

# tstart_analysis = [tstart_1, tstart_2, tstart_3, tstart_4, tstart_5, tstart_6]
# tend_analysis = [tend_1, tend_2, tend_3, tend_4, tend_5, tend_6]

tstart_1 = '2023-04-22'
tend_1 = '2023-05-05'

tstart_2 = '2023-05-06'
tend_2 = '2023-05-20'

tstart_3 = '2023-05-21'
tend_3 = '2023-06-04'

tstart_4 = '2023-06-05'
tend_4 = '2023-06-18'

tstart_5 = '2023-07-04'
tend_5 = '2023-07-18'

tstart_6 = '2023-07-19'
tend_6 = '2023-08-02'

tstart_7 = '2023-08-03'
tend_7 = '2023-08-17'

tstart_8 = '2023-08-18'
tend_8 = '2023-09-01'

tstart_9 = '2023-09-16'
tend_9 = '2023-09-30'

tstart_10 = '2023-10-07'
tend_10 = '2023-10-21'

tstart_11 = '2023-10-22'
tend_11 = '2023-11-05'

tstart_12 = '2023-11-06'
tend_12 = '2023-11-20'

tstart_13 = '2023-11-21'
tend_13 = '2023-12-05'

tstart_all = tstart_1
tend_all = tend_13
tstart_analysis = [tstart_1, tstart_2, tstart_3, tstart_4, tstart_5, tstart_6,
                    tstart_7, tstart_8, tstart_9, tstart_10, tstart_11, tstart_12,
                    tstart_13, tstart_all]
tend_analysis = [tend_1, tend_2, tend_3, tend_4, tend_5, tend_6, tend_7, tend_8,
                  tend_9, tend_10, tend_11, tend_12, tend_13, tend_all]

tstart_analysis = [tstart_all]
tend_analysis = [tend_all]

if config['ANALYSIS_OPTIONS']['budget_bar_plot'] == 'True':
    print('-----------------------------------------------------------------------------')
    print('Emission and deposition bar plot budget')
    print('-----------------------------------------------------------------------------')
    if not os.path.exists(config['FILES_PATHS']['output_results_folder']+'\\Budget'):
        os.makedirs(config['FILES_PATHS']['output_results_folder']+'\\Budget')
    plt.close('all')
    for start, end in zip(tstart_analysis, tend_analysis):
        budget = budget_bar_plot(ds_cal, 'OH_reactivity', start, end, config,
                                  option_mean=2, option_barplot=2,
                                  y_min=1E-6, y_max=1.5E-1,
                                  save_fig = False)

# =============================================================================
# %% Stacked bar plot
# =============================================================================
# Chose dates
tstart_1 = '2023-04-22'
tend_1 = '2023-05-05'

tstart_2 = '2023-05-06'
tend_2 = '2023-05-20'

tstart_3 = '2023-05-21'
tend_3 = '2023-06-04'

tstart_4 = '2023-06-05'
tend_4 = '2023-06-18'

tstart_5 = '2023-07-04'
tend_5 = '2023-07-18'

tstart_6 = '2023-07-19'
tend_6 = '2023-08-02'

tstart_7 = '2023-08-03'
tend_7 = '2023-08-17'

tstart_8 = '2023-08-18'
tend_8 = '2023-09-01'

tstart_9 = '2023-09-16'
tend_9 = '2023-09-30'

tstart_10 = '2023-10-07'
tend_10 = '2023-10-21'

tstart_11 = '2023-10-22'
tend_11 = '2023-11-05'

tstart_12 = '2023-11-06'
tend_12 = '2023-11-20'

tstart_13 = '2023-11-21'
tend_13 = '2023-12-05'

tstart_analysis = [tstart_1, tstart_2, tstart_3, tstart_4, tstart_5, tstart_6,
                    tstart_7, tstart_8, tstart_9, tstart_10, tstart_11, tstart_12]
tend_analysis = [tend_1, tend_2, tend_3, tend_4, tend_5, tend_6, tend_7, tend_8,
                  tend_9, tend_10, tend_11, tend_12]

# tstart_analysis = [tstart_1, tstart_2, tstart_4, tstart_9, tstart_10, tstart_11]
# tend_analysis = [tend_1, tend_2, tend_4, tend_9, tend_10, tend_11]

# mz_list to plot
# mz_list_flux = [69.07, 137.133, 33.034, 135.117, 47.049, 47.013, 61.028, 71.049, 87.08, 75.044,
#            45.034, 31.018, 83.049, 101.06, 139.112, 153.055, 205.195, 93.07, 153.127]

# mz_list_flux_bis = [69.07, 137.13, 33.034, 135.12, 47.049, 47.013, 57.07, 59.049, 85.065, 99.08,
#                     61.028, 89.06, 45.034, 83.086, 87.08, 111.08, 31.018 ,71.049, 75.044 ,73.065, 101.06,
#                     83.049, 139.11, 205.2, 43.054, 153.05, 153.13, 69.034, 151.06, 93.07, 121.1, 107.09]

mz_list_flux = [69.07, 137.13 ,33.034, 135.12, 59.049, 47.049, 61.028, 45.034, 83.086, 31.018,
                    47.013, 87.08, 71.049, 73.065, 75.044, 101.06, 57.07, 83.049, 139.11, 153.05, 205.2,
                    111.08, 43.054, 85.065, 89.06, 99.08, 153.13, 93.07, 121.1, 151.06, 69.034, 49.011, 107.09]

mz_list_flux_OH = [69.07, 137.13, 83.086, 87.08, 31.018, 71.049, 73.065, 47.049, 135.12, 45.034,
                   83.049, 101.06, 33.034, 139.11, 205.2, 57.07, 59.049, 85.065, 99.08, 43.054, 69.034,
                   151.06, 47.013, 61.028, 111.08, 75.044, 89.06 ,93.07, 121.1 ,107.09, 153.13, 153.05]

mz_list_dep = [31.018, 33.034, 45.034, 47.013, 47.049, 59.049, 61.028, 63.044]

if config['ANALYSIS_OPTIONS']['budget_stacked_plot'] == 'True':
    print('-----------------------------------------------------------------------------')
    print('Emission and deposition stacked plot budget')
    print('-----------------------------------------------------------------------------')
    if not os.path.exists(config['FILES_PATHS']['output_results_folder']+'\\Budget'):
        os.makedirs(config['FILES_PATHS']['output_results_folder']+'\\Budget')

    list_budget_hour = budget_stacked_plot(ds=ds_cal,
                        var_plot='flux_cal_ug',
                        tstart_list=tstart_analysis,
                        tend_list=tend_analysis,
                        mz_list=mz_list_flux,
                        n_comp=10,
                        only_sign=None,
                        share_y=False,
                        save_fig=True,
                        directory=config['FILES_PATHS']['output_results_folder']+'\\Budget')

# =============================================================================
# %% Stacked bar plot (for poster)
# =============================================================================
# tstart_plot = ['2023-04-22', '2023-05-06', '2023-06-03', '2023-09-23', '2023-10-07', '2023-10-21']
# tend_plot = ['2023-05-06', '2023-05-20', '2023-06-17', '2023-10-07', '2023-10-21', '2023-11-04']

# tstart_plot = ['2023-05-06']
# tend_plot = ['2023-05-20']
# # tstart_plot = [tstart_1]
# # tend_plot = [tend_1]

# mz_list_flux = [69.07, 137.13 ,33.034, 135.12, 59.049, 47.049, 61.028, 45.034, 83.086, 31.018,
#                     47.013, 87.08, 71.049, 73.065, 75.044, 101.06, 57.07, 83.049, 139.11, 153.05, 205.2,
#                     111.08, 43.054, 85.065, 89.06, 99.08, 153.13, 93.07, 121.1, 151.06, 69.034, 49.011, 107.09]

# mz_list_flux_OH = [69.07, 137.13, 83.086, 87.08, 31.018, 71.049, 73.065, 47.049, 135.12, 45.034,
#                     83.049, 101.06, 33.034, 139.11, 205.2, 57.07, 59.049, 85.065, 99.08, 43.054, 69.034,
#                     151.06, 47.013, 61.028, 111.08, 75.044, 89.06 ,93.07, 121.1 ,107.09, 153.13, 153.05]

# mz_list = [mz_list_flux, mz_list_flux_OH]
# var_plot_list = ['flux_cal_ug', 'flux_OH_reactivity']
# unit_list = ['Mean flux [$µg$ $m^{-2}$ $s^{-1}$]', 'Mean OH reactivity flux [$m$ $s^{-2}$]']
# loc_legend_list = ['upper left', 'upper right']
# color_list = [['darkblue', 'royalblue', 'lightskyblue', 'aqua', 'springgreen',
#                 'forestgreen', 'olive', 'gold', 'orange', 'saddlebrown',
#                 'black', 'red'],
#               ['darkblue', 'royalblue', 'orange', 'teal', 'saddlebrown',
#                 'magenta', 'purple', 'forestgreen', 'aqua', 'gold',
#                 'black', 'red']]

# # color_list = [['#f2bfbf', '#d6616b', '#cedb9c', '#8ba252', '#de9ed5',
# #                 '#a55194', '#9c9ede', '#5253a2', '#fdd0a2', '#fd973c',
# #                 'black', 'red'],
# #               ['#f2bfbf', '#d6616b', '#fdd0a2', 'black', '#fd973c',
# #                 'black', 'black', '#a55194', '#8ba252', '#5253a2',
# #                 'black', 'red']]

# n_comp = 10
# cm = 1/2.54

# for tstart, tend in zip(tstart_plot, tend_plot):

#     fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=False,
#                             figsize=(20*cm, 10*cm))
#     plt.style.use('seaborn-whitegrid')
#     fig.patch.set_facecolor('#fffbe9')
    
#     for i, ax, mz_list_plot, var_plot, unit, loc_legend, color \
#         in zip(range(0, 2), axs.ravel(), mz_list, var_plot_list, unit_list, loc_legend_list, color_list):
#         # Select the mz to plot
#         ds_sel = ds_cal.sel(time=slice(tstart, tend)).sel(mz=mz_list_plot, method='nearest')
        
#         # Compute mean flux per hour
#         budget_hour = (ds_sel[var_plot].groupby('time.hour').mean().rename({'hour': 'time'})\
#                     .reindex(time=np.arange(24))).to_pandas()
                   
#         # Round mz values in column names to have 3 digits after decimal point
#         budget_hour.columns = [f"{col:.3f}" if isinstance(col, (float, np.floating)) else col for col in budget_hour.columns]
        
#         # Ensure 5 significant digits for mz values
#         budget_hour.columns = [f"{col:.5g}" if isinstance(col, (float, np.floating)) else col for col in budget_hour.columns]
        
#         # Only keep the n selected most exchanged and compute the sum of positive
#         # and negative fluxes for the others
#         def sum_pos_neg(row):
#             pos_sum = row.iloc[n_comp:][row.iloc[n_comp:] > 0].sum()
#             neg_sum = row.iloc[n_comp:][row.iloc[n_comp:] < 0].sum()
#             return pd.Series({'other pos': pos_sum, 'other neg': neg_sum})
        
#         budget_hour[['other pos', 'other neg']] = budget_hour.apply(sum_pos_neg, axis=1)
#         budget_hour = budget_hour.iloc[:, :n_comp].join(budget_hour[['other pos', 'other neg']])
        
#         # Plot
#         budget_hour.plot(kind='bar', stacked=True, width=0.9, ax=ax, legend=False,
#                           color=color)
    
#         # Formatting
#         # Chose subplot background color
#         ax.set_facecolor('#fffbe9')
#         # y axis
#         ax.set_ylabel(unit, fontsize=11)
#         ax.tick_params(axis='y', labelsize=10)
#         if i==1:
#             ax.yaxis.tick_right()
#             ax.yaxis.set_label_position('right')
        
#         # x axis
#         ax.set_xticks([0, 6, 12, 18])
#         ax.set_xticklabels(['00:00', '06:00', '12:00', '18:00'], rotation=0)
#         ax.tick_params(axis='x', labelsize=11)
#         ax.set_xlabel('', rotation=0)
        
#         # legend, title
#         legend = ax.legend(labels=budget_hour.columns, loc=loc_legend, fontsize=7)
#         plt.subplots_adjust(wspace=0.05)
#         plt.suptitle(datetime.strptime(tstart, '%Y-%m-%d').strftime('%d-%b-%Y') +
#                       ' to ' + datetime.strptime(tend, '%Y-%m-%d').strftime('%d-%b-%Y'),
#                       fontsize=11, y=0.92)
        
#         # Epurate
#         ax.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5)
#         ax.set_axisbelow(True)
    
#         os.chdir(r'C:\Users\u230638\DOC2\2 Recherche\5 Communication\1 Présentations orales\2024 06 GRS-GRC\Poster\Images')
#         plt.savefig('BVOC_budget_' + tstart + '_' + tend + '.png', transparent=False)
#         # plt.close('all')

# =============================================================================
# Create figure legend
# =============================================================================

# Lists
# color_list = ['darkblue', 'royalblue', 'lightskyblue', 'aqua', 'springgreen', 'forestgreen',
#               'olive', 'gold', 'orange', 'saddlebrown', 'teal', 'magenta', 'purple']
# mz_list_flux = ['69.070', '137.133', '33.034', '135.117', '59.049', '47.049', '61.028', '45.034', '83.086', '31.018', '87.080', '71.049', '73.065']
# name_list = ['isoprene', 'MTs', 'methanol', 'p-cymene', 'acetone', 'ethanol', 'acetic acid',
#               'acetaldehyde', 'GLV fragment', 'formaldehyde', 'MBO', 'MACR+MVK', 'MEK']

# color_list = ['darkblue', 'royalblue', 'lightskyblue', 'aqua', 'springgreen',
#               'forestgreen', 'olive', 'gold', 'orange', 'saddlebrown',
#               'black', 'red']
# mz_list_flux = ['69.070', '137.133', '33.034', '135.117', '59.049', '47.049',
#                 '61.028', '45.034', '83.086', '31.018', 'other > 0', 'other < 0']
# name_list = ['isoprene', 'MTs', 'methanol', 'p-cymene', 'acetone', 'ethanol', 'acetic acid',
#               'acetaldehyde', 'GLV fragment', 'formaldehyde', '', '']

# # Function to add space between columns
# def format_label(mz, name, space=1):
#     return f"{mz:<{max_mz_length}}{' ' * space}{name}"

# # Create fig
# fig, ax = plt.subplots(figsize=(2.5, len(name_list) * 0.25))
# # fig.patch.set_facecolor('#fffbe9')
# ax.axis('off')

# # Compute max length
# max_mz_length = max(len(mz) for mz in mz_list_flux)  # Longueur max des valeurs mz

# # Patch with expected space
# space_between = 3  # Ajustez ce nombre pour augmenter l'espacement
# patches = [mpatches.Patch(color=color_list[i], label=format_label(mz_list_flux[i], name_list[i], space_between)) for i in range(len(name_list))]

# # Legend
# ax.legend(handles=patches, loc='center', frameon=False)

# # Remove space
# plt.subplots_adjust(left=0., right=1., top=1., bottom=0.)
# plt.tight_layout()

# os.chdir(r'C:\Users\u230638\DOC2\2 Recherche\4 Analyse donnees\VOCs\EC\3 Visualisation\Output\2023\Top\All\Budget')
# plt.savefig('legend.png', transparent=False)

# =============================================================================
# %% Pie plot
# =============================================================================
# Chose dates
tstart_1 = '2023-04-22'
tend_1 = '2023-05-05'

tstart_2 = '2023-05-06'
tend_2 = '2023-05-20'

tstart_3 = '2023-05-21'
tend_3 = '2023-06-04'

tstart_4 = '2023-06-05'
tend_4 = '2023-06-18'

tstart_5 = '2023-07-04'
tend_5 = '2023-07-18'

tstart_6 = '2023-07-19'
tend_6 = '2023-08-02'

tstart_7 = '2023-08-03'
tend_7 = '2023-08-17'

tstart_8 = '2023-08-18'
tend_8 = '2023-09-01'

tstart_9 = '2023-09-16'
tend_9 = '2023-09-30'

tstart_10 = '2023-10-07'
tend_10 = '2023-10-21'

tstart_11 = '2023-10-22'
tend_11 = '2023-11-05'

tstart_12 = '2023-11-06'
tend_12 = '2023-11-20'

tstart_13 = '2023-11-21'
tend_13 = '2023-12-05'

tstart_analysis = [tstart_1, tstart_2, tstart_3, tstart_4, tstart_5, tstart_6,
                    tstart_7, tstart_8, tstart_9, tstart_10, tstart_11, tstart_12]
tend_analysis = [tend_1, tend_2, tend_3, tend_4, tend_5, tend_6, tend_7, tend_8,
                  tend_9, tend_10, tend_11, tend_12]

if config['ANALYSIS_OPTIONS']['plot_pie'] == 'True':
    print('-----------------------------------------------------------------------------')
    print('Pie chart')
    print('-----------------------------------------------------------------------------')
    if not os.path.exists(config['FILES_PATHS']['output_results_folder']+'\\Pie_plot'):
        os.makedirs(config['FILES_PATHS']['output_results_folder']+'\\Pie_plot')
    
    list_budget = pie_plot(ds=ds_cal,
                           mz_list=ds_cal.mz.to_pandas().to_list(),
                           var_plot='OH_reactivity',
                           tstart_list=tstart_analysis,
                           tend_list=tend_analysis,
                           range_plot='both',
                           save_fig=True,
                           directory=config['FILES_PATHS']['output_results_folder']+'\\Pie_plot')

# =============================================================================
# %% Diel multiplot over a period
# =============================================================================
#  Diel multiplot evolution only programmed if fluxes calibrated
if config['ANALYSIS_OPTIONS']['apply_calib'] == 'True' and\
    config['ANALYSIS_OPTIONS']['multiplot_diel_period'] == 'True':
        print('-----------------------------------------------------------------------------')
        print('Diel multiplot creation')
        print('-----------------------------------------------------------------------------')
        # Create a new folder to store results
        if not os.path.exists(config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Multiplot'):
            os.makedirs(config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Multiplot')

        # tstart_analysis = [tstart_1, tstart_2, tstart_3, tstart_4, tstart_5, tstart_6]
        # tend_analysis = [tend_1, tend_2, tend_3, tend_4, tend_5, tend_6]
        tstart_analysis = [tstart_1, tstart_2, tstart_3, tstart_4, tstart_5, tstart_6,
                           tstart_7, tstart_8, tstart_9, tstart_10, tstart_11, tstart_12, tstart_13]
        tend_analysis = [tend_1, tend_2, tend_3, tend_4, tend_5, tend_6, tend_7, tend_8,
                         tend_9, tend_10, tend_11, tend_12, tend_13]
        for start, end in zip(tstart_analysis, tend_analysis):
            multiplot_ds_diel_period(ds=ds_cal,
                mz_list='all',
                tstart=start,
                tend=end,
                concat_option=config['USER_OPTIONS']['concat_option'],
                mz_attributed=False,
                suptitle='Diel fluxes',
                var_diel='flux_cal_ug',
                var_right_met=True,
                mz_info='mz_chemform',
                ylim=None,
                hline=None,
                numerical_info=None,
                figclose=True,
                figwidth=min(15, ds.mz.size*4/3),
                figsave=True,
                output_folder=config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Multiplot')

# =============================================================================
# %% Diel plot per period
# =============================================================================
#  Diel plot evolution only programmed if fluxes calibrated
if config['ANALYSIS_OPTIONS']['apply_calib'] == 'True' and\
    config['ANALYSIS_OPTIONS']['plot_diel_period'] == 'True':
        print('-----------------------------------------------------------------------------')
        print('Diel plot creation')
        print('-----------------------------------------------------------------------------')
        # Create a new folder to store results
        if not os.path.exists(config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Plot'):
            os.makedirs(config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Plot')
        list_incr = [14]
        plt.style.use('seaborn-whitegrid')
        # list_mz_all = ds.mz.values.tolist()
        # list_mz_all = [list_mz_all[i:i+4] for i in range(0, len(list_mz_all), 4)]
        list_mz_cal = [[33.03, 69.069, 137.13],
                    [33.03, 47.049, 59.049],
                    [31.018, 45.033, 69.033],
                    [47.013, 61.028],
                    [69.069, 71.049, 75.044, 101.06],
                    [93.0707, 107.087, 121.101],
                    [137.132, 139.112, 153.127, 135.117],
                    [151.112, 153.055, 205.195],
                    [43.05, 73.065],
                    [83.049, 87.08]]  
        
        list_mz_cal = ds_cal.mz.values.tolist()
        list_mz_cal = [list_mz_cal[i:i+4] for i in range(0, len(list_mz_cal), 4)]       

        # list_mz_sign_uncal = [[17.027, 42.034, 44.997, 57.034],
        #                         [63.044, 68.049, 68.062, 77.06],
        #                         [80.049, 82.041, 92.062, 94.065],
        #                         [96.093, 97.101, 101.1, 138.05],
        #                         [143.14, 149.12, 195.07, 203.94],
        #                         [329.84]]       
        # list_mz_uncal = [[29.997, 55.054, 70.041, 71.086],
        #                   [72.057, 73.028, 74.06, 79.054],
        #                   [81.03, 81.045, 84.093, 87.044],
        #                   [90.948, 91.039, 91.054, 94.041],
        #                   [97.065, 99.044, 103.08, 105.05],
        #                   [111.12, 115.08, 115.11, 137.06],
        #                   [140.12, 141.05, 143.11]]
        
        list_mz_cal_bis = [[57.07, 83.086, 85.065],
                            [89.06, 99.08, 111.08]]
        
        list_mz_isop_MT = [[43.054, 57.034, 57.07, 68.049, 69.07, 71.049, 71.086,
                            75.044, 77.06, 80.049, 81.07, 83.049, 85.065, 87.044,
                            87.08, 92.062, 93.07, 94.065, 96.093, 101.06, 111.08,
                            115.08, 135.12, 137.13, 139.11, 140.12, 141.05, 153.05, 153.13, 205.2]]     
        list_mz_methanol = [[31.018, 33.034, 45.034, 47.013, 47.049, 59.049, 61.028, 63.044]]
        list_mz_deposition = [[17.027, 29.997, 44.997, 70.041, 195.07, 329.84]]
        list_mz_other = [[42.034, 81.03, 81.045, 82.041, 137.06]]
        list_mz_budbreak = [[83.086, 55.054, 84.093, 97.101, 101.1, 143.14]]
        list_mz_noisy = [[69.034, 72.057, 73.028, 73.065, 74.06, 79.054, 89.06,
                          90.948, 91.039, 91.054 ,94.041, 97.065, 99.044, 99.08,
                          103.08, 105.05, 111.12, 115.11, 121.1, 138.05, 143.11, 149.12]]
        
        list_mz_poster = [[33.033, 69.069, 83.086]]
               
        list_dates_poster = [['2023-04-22', '2023-05-05'],
                             ['2023-05-07', '2023-05-20'],
                             ['2023-06-03', '2023-06-17'],
                             ['2023-09-23', '2023-10-06'],
                             ['2023-10-08', '2023-10-20'],
                             ['2023-10-22', '2023-11-04']]

        # list_color_poster = ['royalblue', 'aqua', 'springgreen', 'forestgreen', 'olive', 'orange']
        
        list_color_poster = ['#1d4458', '#952f4a', '#1d4458', '#952f4a', '#1d4458', '#952f4a']
        
        # all_mz_lists = [list_mz_budbreak, list_mz_methanol, list_mz_deposition, list_mz_other]
        # all_mz_names = ['list_mz_budbreak', 'list_mz_methanol', 'list_mz_deposition', 'list_mz_other']
        # for list_mz, folder_subname in zip(all_mz_lists, all_mz_names):
        #     if not os.path.exists(config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Plot'+'\\'+folder_subname):
        #         os.makedirs(config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Plot'+'\\'+folder_subname)
        #     for i in list_incr:
        #         plot_ds_diel_period(ds, config,
        #                               tstart1='2023-04-24',
        #                               tlast='2023-11-30',
        #                               bud_break='2023-05-01',
        #                               leaf_fall='2023-10-31',
        #                               var_diel=['flux', 'TA_1_1_1', 'RH_1_1_1', 'WD_1_1_1'],
        #                               var_diel_std=['flux_std', 'TA_1_1_1_std', 'RH_1_1_1_std', 'WD_1_1_1_std'],
        #                               single_mean=True,
        #                               mz_list=list_mz,
        #                               incr=i,
        #                               figsave=True,
        #                               output_folder=config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Plot'+'\\'+folder_subname)
        for i in list_incr:
            # plot_ds_diel_period(ds_cal, config,
            #                     tstart1='2023-04-22',
            #                     tlast='2023-11-30',
            #                     bud_break='2023-05-01',
            #                     leaf_fall='2023-10-31',
            #                     var_diel=['flux_cal_ug', 'TA_1_1_1', 'RH_1_1_1'],
            #                     var_diel_std=['flux_cal_ug_std', 'TA_1_1_1_std', 'RH_1_1_1_std'],
            #                     single_mean=True,
            #                     mz_list=list_mz_poster,
            #                     incr=i,
            #                     figsave=True,
            #                     output_folder=config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Plot')
            plot_ds_diel_period(ds_cal, config,
                                tstart1='2023-04-22',
                                tlast='2023-11-30',
                                bud_break=None,
                                leaf_fall=None,
                                var_diel=['flux_cal_ug', 'TA_1_1_1', 'RH_1_1_1'],
                                var_diel_std=['flux_cal_ug_std', 'TA_1_1_1_std', 'RH_1_1_1_std'],
                                list_label=['$T_{air}$ [°C]', 'RH [%]',
                                            'Methanol flux\n[$µg$ $m^{-2}$ $s^{-1}$]',
                                            'Isoprene flux\n[$µg$ $m^{-2}$ $s^{-1}$]',
                                            'GLV fragment flux\n[$µg$ $m^{-2}$ $s^{-1}$]'],
                                list_date=None,
                                list_color=list_color_poster,
                                single_mean=True,
                                mz_list=list_mz_poster,
                                incr=i,
                                title=False,
                                figsave=True,
                                output_folder=config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Plot')



# plot_ds_diel_period(ds, config,
#                     tstart1='2023-04-24',
#                     tlast='2023-11-30',
#                     bud_break='2023-05-01',
#                     leaf_fall='2023-10-31',
#                     var_diel=['dep_velocity', 'TA_1_1_1', 'RH_1_1_1'],
#                     var_diel_std=['dep_velocity_std', 'TA_1_1_1_std', 'RH_1_1_1_std'],
#                     single_mean=True,
#                     mz_list=[[33.03]],
#                     incr=7,
#                     figsave=False,
#                     output_folder=config['FILES_PATHS']['output_results_folder']+'\\Diel_plots'+'\\Plot')

# =============================================================================
# %% Clock plot
# =============================================================================
# print('-----------------------------------------------------------------------------')
# print('Clock plot (variable as a function of hour)')
# print('-----------------------------------------------------------------------------')
# plot_ds_clock(ds=ds_cal,
#               var='flux_cal_ug',
#               mz_stat=47.013,
#               min_max='idxmin',
#               freq='D')
#

# =============================================================================
# %% Relative contribution to the total flux
# =============================================================================
# print('-----------------------------------------------------------------------------')
# print('Relative contribution to the total flux')
# print('-----------------------------------------------------------------------------')
# # Options
# mz_list_sign = [31.018, 33.034, 42.034, 43.054, 45.034, 47.013, 47.049, 57.07,
#                 59.049, 61.028, 69.034, 69.07, 71.049, 73.065, 75.044, 83.049,
#                 85.065, 87.08, 89.06, 99.08, 101.06, 111.08, 121.1, 135.12, 137.13,
#                 139.11, 153.05, 153.13, 205.2, 63.044, 70.041, 72.057, 73.028,
#                 74.06, 77.06, 79.054, 80.049, 81.045, 82.041, 87.044, 91.039, 92.062,
#                 94.041, 94.065, 97.065, 99.044, 101.1, 103.08, 105.05, 115.08, 115.11,
#                 137.06, 138.05, 140.12, 141.05, 143.11, 143.14, 149.12]

# tstart_analysis = ['2023-04-22']
# tend_analysis = ['2023-12-05']

# relative_contribution(ds=ds_cal_kinetic,
#                       var_plot='flux_cal_ug',
#                       mz_list=mz_list_sign,
#                       tstart_list=tstart_analysis,
#                       tend_list=tend_analysis)

# =============================================================================
# %% Mass defect plot
# =============================================================================
# print('-----------------------------------------------------------------------------')
# print('Mass defect plot')
# print('-----------------------------------------------------------------------------')
# kendrick_plot(ds_cal_kinetic,
#               plot_type='KMD',
#               scale_conc=True,
#               annotate_all=False,
#               show_lines=True)

# =============================================================================
# %% BVOC fingerprint
# =============================================================================
# print('-----------------------------------------------------------------------------')
# print('BVOC fingerprint concentration = f(mz)')
# print('-----------------------------------------------------------------------------')
# bvoc_fingerprint(ds)

# =============================================================================
# %% Ratio of mz1 to mz2 fluxes
# =============================================================================
# print('-----------------------------------------------------------------------------')
# print('Ratio of mz1 to mz2 fluxes')
# print('-----------------------------------------------------------------------------')
# daily_value = 'max'  # 'mean' or 'max'

# # Compute ds_mean or ds_max
# if daily_value == 'mean':
#     if 'ds_mean' not in globals():
#         ds_mean = ds.resample(time='D').mean()
#     ds_plot = ds_mean
# elif daily_value == 'max':
#     if 'ds_max' not in globals():
#         ds_max = ds.resample(time='D').max()
#     ds_plot = ds_max

# plot_ratio_mz(ds = ds_plot,
#               mz1 = 33.033,
#               mz2 = 69.069,
#               plot_type = daily_value)