# -*- coding: utf-8 -*-
"""
Created on Tue May 23 10:38:09 2023

@author: u230638
"""

import pandas as pd
import numpy as np
import os
import xarray as xr
import math
from TOFEC_concat_hdf5_IRGA import concat_hdf5_old, concat_hdf5_new
from TOFEC_plot_xarray import plot_ds_time

# =============================================================================
# %% Extract the attributes for a given mz
# =============================================================================


def mz_attr(ds, mz_info, mz_interest):
    """
    Used to extract from the xarray Dataset the different IDs of a given mz

    Parameters
    ----------
    ds : xarray Dataset from which the info must be extracted.
    mz_info : information required concerning a given mz
    mz_interest : for which mz the information is asked. Set to 'all' or a given mz.

    Returns
    -------
    attrs : array of different IDs for each mz.

    """
    if mz_interest == 'all':
        attrs = pd.DataFrame([row.split("/") for row in ds[mz_info].to_numpy()])
        attrs.index = np.around(ds.mz.values, 2)
    else:
        attrs = pd.DataFrame(ds[mz_info].isel(mz=np.abs(ds.mz - mz_interest).argmin())
                             .to_numpy()[()].split("/"))

    return attrs

# =============================================================================
# %% Verify the concatenation of GEddySoft daily output files and start it if necessary
# =============================================================================


def verif_concat(config):
    """
    Verify if the concatenation of GEddySoft daily output files has already been performed.
    If not, launch TOFEC_concat_hdf5.py

    Parameters
    ----------
    config : file containing files paths, files names and options

    Returns
    -------
    None.

    """

    # Verifiy if hdf5 result files have already been concatenated
    if config['USER_OPTIONS']['output_type'] == 'old':       
        hdf5_name = 'results_' + eval(config['USER_OPTIONS']['date_files_selection'])[0][0:10]\
            + '_' + eval(config['USER_OPTIONS']['date_files_selection'])[1][5:10]\
                + '_' + config['USER_OPTIONS']['concat_option'] + '.nc'
    else:
        hdf5_name = config['FILES_NAMES']['GEddySoft_results_files_prefix']\
            + eval(config['USER_OPTIONS']['date_files_selection'])[0][0:10] + '_' \
                + eval(config['USER_OPTIONS']['date_files_selection'])[1][5:10] + '.nc'
    if os.path.exists(os.path.join(config['FILES_PATHS']['concat_results_output_folder'], hdf5_name)):
        print('-----------------------------------------------------------------------------')
        print('GEddySoft daily hdf5 outputs have already been concatenated for the selected\n' +
              'time period')
        print('-----------------------------------------------------------------------------')
    # If not, launch TOFEC_concat_hdf5.py to concatenate daily results files
    else:
        print('-----------------------------------------------------------------------------')
        print('GEddySoft daily hdf5 outputs have to be concatenated for the selected\n' +
              'time period')
        print('The TOFEC_concat_hdf5 routine will start now. This task may take a few ' +
              'minutes.\nA message will be displayed once the procedure over.')
        print('-----------------------------------------------------------------------------')
        print('TOFEC_concat_hdf5.py running...')
        
        if config['USER_OPTIONS']['output_type'] == 'old':
            concat_hdf5_old(config)
        else:
            concat_hdf5_new(config)

        print('TOFEC_concat_hdf5.py executed with success.')

# =============================================================================
# %% Compute xh-means
# =============================================================================



def mean_xh(ds, variable, freq, config):
    """
    

    Parameters
    ----------
    ds : input xarray dataset
    var : variable to be rounded
    freq : frequency to be rounded ('2H', '3H', '4H', etc.)
    config : file containing files paths, files names and options

    Returns
    -------
    ds_xh_mean : xarray dataset with variables computed to xh means

    """
    # # Compute 2h means on selected variable except mz information
    # if config['USER_OPTIONS']['concat_option'] == 'GLOVOCS':
    #     ds_xh_mean = xr.merge([ds.drop_vars(['mz_glov', 'mz_match', 'mz_name',
    #                                          'mz_chemform']).resample(time='2H').mean(),
    #                            ds.get(['mz_glov', 'mz_match', 'mz_name', 'mz_chemform'])])
    # elif config['USER_OPTIONS']['concat_option'] == 'IONICON':
    #     ds_xh_mean = xr.merge([ds.drop_vars(['mz_chemform']).resample(time='2H').mean(),
    #                            ds.get(['mz_chemform'])])
    print('-----------------------------------------------------------------------------')
    print('Computation of ' + freq + ' mean ds')
    print('-----------------------------------------------------------------------------')
    
    # Remove the time dimension of ds['mz_chemform'] if needed   
    if 'mz_chemform' in ds and ds['mz_chemform'].dims.count('time') > 0 and ds['mz_chemform'].sizes['time'] > 1:
        ds['mz_chemform'] = ds['mz_chemform'].isel(time=0).squeeze()
        
    # Compute 2h means on selected variable except mz information
    if config['USER_OPTIONS']['output_type'] == 'old':
        if config['USER_OPTIONS']['concat_option'] == 'GLOVOCS':
            ds_xh_mean = xr.merge([xr.merge([ds[variable], ds['calcoef_norm']]).resample(time=freq, skipna=True).mean().compute(),
                                   ds.get(['mz_glov', 'mz_match', 'mz_name', 'mz_chemform'])])
        elif config['USER_OPTIONS']['concat_option'] == 'IONICON':
            ds_xh_mean = xr.merge([xr.merge([ds[variable], ds['calcoef_norm']]).resample(time=freq, skipna=True).mean().compute(),
                                   ds.get(['mz_chemform'])])
    else:
            ds_xh_mean = ds[variable].resample(time=freq, skipna=True).mean().compute()

    # Keep the attributes in the grouped dataset
    for var in ds_xh_mean:
        ds_xh_mean[var].attrs = ds[var].attrs
    return ds_xh_mean

# =============================================================================
# %% Create multiplots
# =============================================================================


def multiplot(ds, config, output_folder):
    """


    Parameters
    ----------
    ds : input xarray dataset
    config : file containing files paths, files names and options
    output_folder : folder to store all the multiplot figures

    Returns
    -------
    Multiplot of selected variables for all m/z values.

    """
    
    for i in eval(config['ANALYSIS_OPTIONS']['multiplot_var']):
        variable = [i]
        print('variable plotted: ', variable)
        
        plot_type = eval(config['ANALYSIS_OPTIONS']['multiplot_var'])[i]
        
        ds_plot = mean_xh(ds, variable, '2H', config)
        
        
        plot_ds_time(ds=ds_plot,
            mz_list='all',
            tstart=None,
            tend=None,
            concat_option=None,
            mz_attributed=False,
            suptitle=i,
            var_left=variable,
            kind_left=plot_type,
            var_right=None,
            kind_right='scatter',
            var_right_met=True,
            mz_info='mz_chemform',
            ylim=None,
            hline=None,
            vline=None,
            numerical_info=None,
            figclose=True,
            figwidth=30,
            figsave=True,
            output_folder=output_folder,
            config=config)