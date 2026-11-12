# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 10:42:53 2022

@author: u230638
"""

import os
import pandas as pd
from datetime import datetime, timedelta
import xarray as xr
import numpy as np
from TOFEC_list_files import list_files_datestring
import sys
import time
import h5py
import re


def concat_hdf5_new(config):
    """
    Concatenate GEddySoft new daily output files into xarray datasets for the desired
    variables. Save it under a .nc data format.
    Parameters
    ----------
    config : file containing files paths, files names and options

    Returns
    -------
    None.

    """

    # =============================================================================
    # List files
    # =============================================================================
    files = list_files_datestring(config['FILES_PATHS']['GEddySoft_results_output_folder'],
                                  config['FILES_NAMES']['GEddySoft_results_files_prefix'],
                                  '.hdf5', 'yyyymmdd',
                                  eval(config['USER_OPTIONS']['date_files_selection']))
    
    # Compute the minimum file size
    files['size'] = [os.path.getsize(os.path.join(path, name))/1000. for path, name in zip(files['path'], files['name'])]
    min_file_size = pd.Series(files['size']).quantile(0.05)
    
    # =============================================================================
    # Prepare file reading
    # =============================================================================
    # List to store all the datasets
    datasets = []
    
    # Loop through all the files in 'files'
    for i in range(len(files['path'])):
        if files['size'][i] < min_file_size:
            print('Not read file: ' + files['name'][i] + ' - File too short')
            continue
        print('Read file: ' + files['name'][i])
        filepath = files['path'][i] + '\\' + files['name'][i]
        
        with h5py.File(filepath, mode='r') as results:
    
            # Get the time and freq values
            time = results['time'][()]
            freq = results['freq'][()]  # New dimension
    
            # Initialize lists to store data
            mz_list = []
            data_vars = {}
            
    # =============================================================================
    #       Read tracer variables
    # =============================================================================
            # Get numeric keys in the TRACER group (e.g., '0', '1', etc.)
            tracer_keys = [k for k in results['IRGA'].keys() if k.isdigit()]
            
            # Variables
            varlist_tracer = eval(config['USER_OPTIONS']['varlist_tracer'])
            varlist_tracer_freq = ['cospec', 'cospec_scaled', 'spec_scaled', 'spec']
            varlist_tracer_time_only = ['cf_lpf', 'default_CC_kinetic']  
            varlist_tracer_qaqc = eval(config['USER_OPTIONS']['varlist_tracer_qaqc'])
    
            # Loop through each numeric tracer group
            for key in tracer_keys:
                group = results['IRGA'][key]
                mz = group['name'][()]
                if mz == b'H2O dry': mz = b'18'
                if mz == b'N2O dry': mz = b'44'
                mz_list.append(mz)
    
                missing_vars = [v for v in varlist_tracer if v not in varlist_tracer_time_only and v not in group.keys()]
                if missing_vars:
                    print(f"[IRGA/{key}] Missing variables: {missing_vars}")
                
                # Tracer variables
                for var in varlist_tracer:
                    if var not in group.keys():
                        continue
                    values = group[var][()]                
                    # Extract attributes 'description' and 'units'
                    description = group[var].attrs.get('description', '')
                    units = group[var].attrs.get('units', '')
                    
                    # Frequency dependent variables (spectrum, co-spectrum)
                    if var in varlist_tracer_freq:
                        # Shape (time, freq) → will become (mz, time, freq)
                        if var not in data_vars:
                            data_vars[var] = {'data': [], 'attrs': {'name': description, 'units': units}}
                        data_vars[var]['data'].append(values)
                        
                    # Frequency undependent variables
                    else:
                        # Shape (time,) → will become (mz, time)
                        if var not in data_vars:
                            data_vars[var] = {'data': [], 'attrs': {'name': description, 'units': units}}
                        data_vars[var]['data'].append(values)
                
                # QAQC variables
                for var in varlist_tracer_qaqc:
                    if var not in group['qaqc'].keys():
                        continue
    
                    # IPT variable: needs to be separated
                    if var.startswith('IPT'):     
                        # Find IPT names in the description
                        description = group['qaqc'][var].attrs.get('description', '')
                        match = re.search(r"\((.*?)\)", description)
                        # If list_IPT not already created:
                        if 'list_IPT' not in locals():
                            if match:
                                list_IPT = [item.strip() for item in match.group(1).split(',')]
                            else:
                                list_IPT = []
                        # Loop to read each IPT separately
                        for i, var_IPT in enumerate(list_IPT):
                            values = group['qaqc'][var][:, i]
                            var_name = var + '_' + var_IPT
                            # print(var_name)
                            description = var + '_' + var_IPT
                            units = '%'
                            if var_name not in data_vars:
                                data_vars[var_name] = {'data': [], 'attrs': {'name': description, 'units': units}}
                            # Add to data_vars (no need to expand since it only depends on time)
                            data_vars[var_name]['data'].append(values)
                    
                    # Other QAQC variables than IPT
                    else:
                        values = group['qaqc'][var][()]
                        # Extract attributes 'description' and 'units'
                        description = group['qaqc'][var].attrs.get('description', '')
                        units = group['qaqc'][var].attrs.get('units', '')  
                        if var not in data_vars:
                            data_vars[var] = {'data': [], 'attrs': {'name': description, 'units': units}}
                        # Add to data_vars (no need to expand since it only depends on time)
                        data_vars[var]['data'].append(values)
                        
            # Convert mz to array
            mz_array = np.array(mz_list)
            
            # Variables in tracer that are not time dependent
            for var in varlist_tracer_time_only:
                if (var not in results['IRGA'].keys()) or (var not in varlist_tracer):
                    continue
                values = results['IRGA'][var][()]
                # Extract attributes 'description' and 'units'
                description = results['IRGA'][var].attrs.get('description', '')
                units = results['IRGA'][var].attrs.get('units', '')
                if var not in data_vars:
                    data_vars[var] = {'data': [], 'attrs': {'name': description, 'units': units}}
                # Add to data_vars (no need to expand since it only depends on time)
                data_vars[var]['data'].append(values)
     
    # =============================================================================
    #       Read met variables
    # =============================================================================
            varlist_met = eval(config['USER_OPTIONS']['varlist_met'])
            varlist_met_freq = ['cospec_wT', 'cospec_wT_scaled', 'spec_T', 'spect_T_scaled']
            varlist_met_qaqc = eval(config['USER_OPTIONS']['varlist_met_qaqc'])
            
            # Regular met variables
            for var in varlist_met:
                if var not in results['MET'].keys():
                    continue
                values = results['MET'][var][()]
                # Extract attributes 'description' and 'units'
                description = results['MET'][var].attrs.get('description', '')
                units = results['MET'][var].attrs.get('units', '')
                if var not in data_vars:
                    data_vars[var] = {'data': [], 'attrs': {'name': description, 'units': units}}
                # Add to data_vars (no need to expand since it only depends on time)
                data_vars[var]['data'].append(values)
            
            # QAQC variables
            for var in varlist_met_qaqc:
                if var not in results['MET']['qaqc'].keys():
                    continue
                
                # IPT variables: need to be separated
                if var.startswith('IPT'):
                    # Find IPT names in the description
                    description = results['MET']['qaqc'][var].attrs.get('description', '')
                    match = re.search(r"\((.*?)\)", description)
                    if match:
                        list_IPT = [item.strip() for item in match.group(1).split(',')]
                    else:
                        list_IPT = []
                    # Loop to read each IPT separately
                    for i, var_IPT in enumerate(list_IPT):
                        values = results['MET']['qaqc'][var][:, i]
                        var_name = var + '_' + var_IPT
                        description = var + '_' + var_IPT
                        units = '%'
                        if var not in data_vars:
                            data_vars[var_name] = {'data': [], 'attrs': {'name': description, 'units': units}}
                        # Add to data_vars (no need to expand since it only depends on time)
                        data_vars[var_name]['data'].append(values)
    
                # Other QAQC variables than IPT
                else:
                    values = results['MET']['qaqc'][var][()]
                    # Extract attributes 'description' and 'units'
                    description = results['MET']['qaqc'][var].attrs.get('description', '')
                    units = results['MET']['qaqc'][var].attrs.get('units', '')
                    if var not in data_vars:
                        data_vars[var] = {'data': [], 'attrs': {'name': description, 'units': units}}
                    # Add to data_vars (no need to expand since it only depends on time)
                    data_vars[var]['data'].append(values)
    
    # =============================================================================
    #       Prepare for xarray
    # =============================================================================
            # Format variables for xarray
            data_vars_xr = {}
            for var, data_dict in data_vars.items():
                # print(var)
                data_np = np.squeeze(np.stack(data_dict['data'], axis=0)) # Flatten data
                attrs = data_dict['attrs']
            
                if var in varlist_tracer_freq:
                    # (mz, time, freq)
                    data_vars_xr[var] = (("time", "mz", "freq"), data_np.transpose(1, 0, 2), attrs)
                elif var in varlist_met_freq:
                    # (time, freq)
                    data_vars_xr[var] = (("time", "freq"), data_np, attrs)
                elif (var in varlist_tracer_time_only) or (var in varlist_met) or (var in varlist_met_qaqc)\
                    or (var.startswith('IPT_T')) or (var.startswith('IPT_w')):
                    # (time)
                    data_vars_xr[var] = (("time",), data_np, attrs)
                elif (var in varlist_tracer) or (var in varlist_tracer_qaqc) or (var.startswith('IPT')):
                    # (mz, time)
                    data_vars_xr[var] = (("time", "mz"), data_np.T, attrs)
                    
            # Decode 'time' from bytes to datetime64 if needed
            if np.issubdtype(time.dtype, np.bytes_):
                time = [t.decode() for t in time]
            elif np.issubdtype(time.dtype, np.str_):
                time = list(time)
            else:
                time = time.astype(str)
            
            # Fix the format 'YYYY-MM-DD HH-MM-SS' → 'YYYY-MM-DD HH:MM:SS'
            time_cleaned = [t[:10] + ' ' + t[11:].replace('-', ':') for t in time]
            
            # Convert to datetime64
            time = np.array(time_cleaned, dtype='datetime64[s]')
            
            # Decode 'mz' from bytes to float if needed
            if mz_array.dtype.kind in {'S', 'O'}:  # string or object
                mz_array = np.array([float(m.decode()) if isinstance(m, bytes) else float(m) for m in mz_array])
    
            # Create the xarray Dataset for this file
            ds_temp = xr.Dataset(
                data_vars_xr,
                coords={"time": time, "mz": mz_array, "freq": freq}
            )
    
            # Append the Dataset to the list
            datasets.append(ds_temp)
            
    # Concatenate all datasets along the 'time' dimension
    ds = xr.concat(datasets, dim='time')
    
    # Sort mz by increasing order
    ds = ds.sortby('mz')
    
    # =============================================================================
    # Read cov files
    # =============================================================================
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        # List files
        files_cov = list_files_datestring(config['FILES_PATHS']['GEddySoft_cov_output_folder'],
                                     config['FILES_NAMES']['GEddySoft_cov_files_prefix'],
                                     '.hdf5', 'yyyymmdd',
                                     eval(config['USER_OPTIONS']['date_files_selection']))
        print(files_cov)
        # Minimum size file
        files_cov['size'] = [os.path.getsize(os.path.join(path, name))/1000. for path, name in zip(files_cov['path'], files_cov['name'])]
        min_file_size_cov = pd.Series(files_cov['size']).quantile(0.05)
        
        # Prepare reading
        # List to store all the datasets
        datasets_cov = []
        
        # Loop through all the files in 'files_cov'
        for i in range(len(files_cov['path'])):
            if files_cov['size'][i] < min_file_size_cov:
                print('Not read cov file: ' + files_cov['name'][i] + ' - File too short')
                continue
            print('Read cov file: ' + files_cov['name'][i])
            filepath_cov = files_cov['path'][i] + '\\' + files_cov['name'][i]
            filepath_tracer = files['path'][i] + '\\' + files['name'][i]
            
            # Open both cov and tracer files (the second one for time information)
            with h5py.File(filepath_cov, 'r') as results_cov, h5py.File(filepath_tracer, 'r') as results_tracer:
        
                # Get the time and freq values
                time = results_tracer['time'][()]
                n_lag = results_cov['IRGA']['0']['cov'][()].shape[1]
                lag = pd.Series(range(0, n_lag))
        
                # Initialize lists to store data
                mz_list = []
                data_vars_cov = {}
                
                # Read
                # Get numeric keys in the TRACER group (e.g., '0', '1', etc.)
                tracer_keys = [k for k in results_cov['IRGA'].keys() if k.isdigit()]
                
                for key in tracer_keys:
                    group = results_cov['IRGA'][key]
                    
                    values = group['cov'][()]
                    # Extract attributes 'description' and 'units'
                    description = group['cov'].attrs.get('description', '')
                    units = group['cov'].attrs.get('units', '')
                    if 'cov' not in data_vars_cov:
                        data_vars_cov['cov'] = {'data': [], 'attrs': {'name': description, 'units': units}}
                    data_vars_cov['cov']['data'].append(values)
            
            # Format variables for xarray
            data_vars_xr_cov = {}
            for var, data_dict in data_vars_cov.items():
                # print(var)
                data_np = np.squeeze(np.stack(data_dict['data'], axis=0)) # Flatten data
                attrs = data_dict['attrs']
                data_vars_xr_cov[var] = (("time", "mz", "lag"), data_np.transpose(1, 0, 2), attrs)
            
            # Decode 'time' from bytes to datetime64 if needed
            if np.issubdtype(time.dtype, np.bytes_):
                time = [t.decode() for t in time]
            elif np.issubdtype(time.dtype, np.str_):
                time = list(time)
            else:
                time = time.astype(str)
            
            # Fix the format 'YYYY-MM-DD HH-MM-SS' → 'YYYY-MM-DD HH:MM:SS'
            time_cleaned = [t[:10] + ' ' + t[11:].replace('-', ':') for t in time]
            
            # Convert to datetime64
            time = np.array(time_cleaned, dtype='datetime64[s]')
            
            # Decode 'mz' from bytes to float if needed
            if mz_array.dtype.kind in {'S', 'O'}:  # string or object
                mz_array = np.array([float(m.decode()) if isinstance(m, bytes) else float(m) for m in mz_array])
        
            # Create the xarray Dataset for this file
            ds_temp_cov = xr.Dataset(
                data_vars_xr_cov,
                coords={"time": time, "mz": mz_array, "lag": lag}
            )
            
            # Append the Dataset to the list
            datasets_cov.append(ds_temp_cov)
            
        # Concatenate all datasets along the 'time' dimension
        ds_cov = xr.concat(datasets_cov, dim='time')
        
        # Sort mz by increasing order
        ds_cov = ds_cov.sortby('mz')
            
    # =============================================================================
    # Save xarray dataset into netcdf
    # =============================================================================
    # Save the xarray in netCDF format
    os.chdir(config['FILES_PATHS']['concat_results_output_folder'])
    
    file_name = config['FILES_NAMES']['GEddySoft_results_files_prefix']\
        + eval(config['USER_OPTIONS']['date_files_selection'])[0][0:10] + '_' \
            + eval(config['USER_OPTIONS']['date_files_selection'])[1][5:10] + '.nc'
    
    # ds.to_netcdf(file_name)
    ds.to_netcdf(file_name, engine="scipy")
    
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        # Save the cov xarray in netCDF format
        os.chdir(config['FILES_PATHS']['geddysoft_cov_output_folder'])
        
        file_name = config['FILES_NAMES']['geddysoft_cov_files_prefix']\
            + eval(config['USER_OPTIONS']['date_files_selection'])[0][0:10] + '_' \
                + eval(config['USER_OPTIONS']['date_files_selection'])[1][5:10] + '.nc'
        
        ds_cov.to_netcdf(file_name)
    
    # =============================================================================
    # Save xarray dataset into csv
    # =============================================================================
    # Output csv if requested
    if config['USER_OPTIONS']['output_csv'] == 'True':
        csv_folder = config['FILES_PATHS']['concat_results_output_folder']\
            + '\\csv_'+ eval(config['USER_OPTIONS']['date_files_selection'])[0][0:10]\
                + '_'+eval(config['USER_OPTIONS']['date_files_selection'])[1][5:10]
        # Verify if the csv folder exists
        if not os.path.exists(csv_folder):
            os.makedirs(csv_folder)
        os.chdir(csv_folder)
        data = {}
        
        # Save met variables : everything except the ones depending on frequency and IPT
        for var in [v for v in varlist_met if v not in varlist_met_freq and not v.startswith('IPT')]:
            data[var + ' ' + str(ds[var].attrs)]=ds[var].values        
        df = pd.DataFrame(data, index=ds['time'].values)
        df.index.name = 'TIMESTAMP_END'
        df.to_csv(config['USER_OPTIONS']['meas_level'] + '_MET.csv', index=True)
        
        # Save tracer variables
        col_title = ['mz_' + element for element in ds['mz'].values.astype(str)]
        for var in [v for v in varlist_tracer if v not in varlist_tracer_freq\
                    and v not in varlist_tracer_time_only]:
            col_title_var = [element + '_' + str(ds[var].attrs)
                              for element in col_title]
            df = pd.DataFrame(ds[var].values, index=ds['time'].values,
                              columns=col_title_var)
            df.index.name = 'TIMESTAMP_END'
            df.to_csv(config['USER_OPTIONS']['meas_level'] + '_TRACER_' + var +'.csv',
                      index=True)
        for var in [v for v in varlist_tracer_qaqc if not v.startswith('IPT')]:
            col_title_var = [element + '_' + str(ds[var].attrs)
                              for element in col_title]
            df = pd.DataFrame(ds[var].values, index=ds['time'].values,
                              columns=col_title_var)
            df.index.name = 'TIMESTAMP_END'
            df.to_csv(config['USER_OPTIONS']['meas_level'] + '_QAQC_' + var +'.csv',
                      index=True)  




def concat_hdf5_old(config):
    """
    Concatenate GEddySoft old daily output files into xarray datasets for the desired
    variables. Save it under a .nc data format.

    Parameters
    ----------
    config : file containing files paths, files names and options

    Returns
    -------
    None.

    """

    # =============================================================================
    # %% Functions
    # =============================================================================
    
    def TRACER_2_3Dnp(results, sel_vars_tracer, sel_vars_qaqc):
        """
        Convert hdf5 results file into a 3D numpy array.
        Dim0 = timestamp
        Dim1 = variable
        Dim2 = mz
    
        Parameters
        ----------
        results : 3D np.
        sel_vars : variables to transfer into the 3D np.
    
        Returns
        -------
        The 3D numpy array where two variables have been added as first columns
        of dim2 : time and tracer name, replicated on dim 3
    
        """
    
        # Create an empty 3D numpy array
        np_res_d = np.empty((len(results['time']), len(sel_vars_tracer)+len(sel_vars_qaqc)+2,
                             len(results['IRGA'])))
    
        # Fill-in the 3D np
        # Add time in col 0 of dim 2
        for i in range(len(results['time'])):
            # Format index as datetime
            time_object = results['time'][i]
            # If time_object is read as bytes, decode it
            if isinstance(time_object, bytes):
                time_object = time_object.decode('utf-8')
            ts = pd.to_datetime(time_object, format='%Y-%m-%d %H-%M-%S')
            np_res_d[i, 0, :] = (time.mktime(ts.timetuple()))
    
        # Add tracer name in col 1 of dim 2
        for i in range(len(results['IRGA'])):
            np_res_d[:, 1, i] = results['IRGA'][str(i)]['name'][...]
    
        # Loop throughout the tracer variables
        for i in range(len(sel_vars_tracer)):
            datatemp = np.linspace(0, len(results['IRGA']['0'][sel_vars_tracer[i]])-1,
                                   num=len(results['IRGA']['0'][sel_vars_tracer[i]]))
            for j in range(len(results['IRGA'])):  # loop on mz
                datatemp = np.c_[datatemp, results['IRGA'][str(j)][sel_vars_tracer[i]][...]]
    
            np_res_d[:, i+2, :] = datatemp[:, 1:]
    
        # Loop throughout the qaqc variables
        IPT_split = eval(config['USER_OPTIONS']['varlist_IPT'])
        for i in range(len(sel_vars_qaqc)):
            if sel_vars_qaqc[i] not in IPT_split:
                datatemp = np.linspace(0, len(results['IRGA']['0']['qaqc'][sel_vars_qaqc[i]])-1,
                                       num=len(results['IRGA']['0']['qaqc'][sel_vars_qaqc[i]]))
                for j in range(len(results['IRGA'])):  # loop on mz
                    datatemp = np.c_[datatemp, results['IRGA'][str(j)]['qaqc']
                                     [sel_vars_qaqc[i]][...]]
        
            else:
                # Fill in with IPT data that are contained inside the IPT matrix
                datatemp = np.linspace(0, len(results['IRGA']['0']['qaqc']['IPT'][:,IPT_split.index(sel_vars_qaqc[i])])-1,
                            num=len(results['IRGA']['0']['qaqc']['IPT'][:,IPT_split.index(sel_vars_qaqc[i])]))
                for j in range(len(results['IRGA'])):  # loop on mz
                    datatemp = np.c_[datatemp, results['IRGA'][str(j)]['qaqc']['IPT']
                                     [:,IPT_split.index(sel_vars_qaqc[i])][...]]

            np_res_d[:, len(eval(config['USER_OPTIONS']['varlist_tracer']))+i+2, :]\
                = datatemp[:, 1:]
 
        
        
        
        # Some timestamps in the hdf5 files are not in the good format. Might give np_res_d
        # arrays larger than 48.
        # Find the positions of wrong timestamps
        list_del = []
        for timestamp in range(len(np_res_d[:, :, :])):
            # Not a half-hour
            if np_res_d[timestamp, 0, 0] % (30*60) != 0:
                list_del.append(timestamp)

        # Delete the wrong HH in np_res_d
        for i in reversed(list_del):
            np_res_d = np.delete(np_res_d, i, axis=0)
        return np_res_d
    
    
    def COV_2_3Dnp(results, cov, n_lag):
        """
        Convert hdf5 cov file into a 3D numpy array.
        Dim0 = timestamp
        Dim1 = lag
        Dim2 = mz
    
        Parameters
        ----------
        results : hdf5 results file
        conv = hdf5 cov file
    
        Returns
        -------
        The 3D numpy array where two variables have been added as first columns
        of dim2 : time and tracer name, replicated on dim 3
    
        """
    
        # Create an empty 3D numpy array
        # np_res_d = np.empty((len(results['time']), int(n_lag)+2, len(cov['IRGA'])))
        np_res_d = np.empty((len(results['time']), int(n_lag)+2, len(cov['IRGA'])))
    
        # Fill-in the 3D np
        # Add time in col 0 of dim 2
        for i in range(len(results['time'])):
            # Format index as datetime
            time_object = results['time'][i]
            # If time_object is read as bytes, decode it
            if isinstance(time_object, bytes):
                time_object = time_object.decode('utf-8')          
            ts = pd.to_datetime(time_object, format='%Y-%m-%d %H-%M-%S')
            
            np_res_d[i, 0, :] = (time.mktime(ts.timetuple()))
    
        # Add tracer name in col 1 of dim 2
        for i in range(len(cov['IRGA'])):
            np_res_d[:, 1, i] = results['IRGA'][str(i)]['name'][...]
    
        # Loop throughout the mz values
        for i in range(len(cov['IRGA'])):
            datatemp = np.c_[cov['IRGA'][str(i)]['cov']]      
            # np_res_d might have a longer time dimension than datatemp (because np_res_d
            # time is based on results time and HH can be missing in HH files. So reduce the
            # size of np_res_d
            np_res_d = np_res_d[0:len(datatemp), :, :]
            np_res_d[:, 2:, i] = datatemp[:, :]
    
        # Some timestamps in the hdf5 files are not in the good format. Might give np_res_d
        # arrays larger than 48.
        # Find the positions of wrong timestamps
        list_del = []
        for timestamp in range(len(np_res_d[:, :, :])):
            # Not a half-hour
            if np_res_d[timestamp, 0, 0] % (30*60) != 0:
                list_del.append(timestamp)
    
        # Delete the wrong HH in np_res_d
        for i in reversed(list_del):
            np_res_d = np.delete(np_res_d, i, axis=0)
  
        return np_res_d
    
    
    def MET_2_2Dnp(results, sel_vars):
        """
        convert the MET part of the results dictionnary to a 2D numpy array
        with dim1=time, dim2=var
        Arguments:
            results: the dictionnary coming from a daily output file of GeddySoft-TOF
            sel_vars: list of strings with the MET sub-keys of the dictionnary to be converted
        Returns:
            The 2D numpy array where one variable have been added as first row of dim2 : time
        Comment:
            Obviously not the most efficient code but uses the same logic as
            dict_TRACER_2_3Dnp
        """
        # put the daily file values in a 2D numpy array
        np_res_d = np.empty((len(results['time']), len(sel_vars)+1))
    
        # add time in col 0 of dim 2
        for i in range(len(results['time'])):
            # Format index as datetime
            time_object = results['time'][i]
            # If time_object is read as bytes, decode it
            if isinstance(time_object, bytes):
                time_object = time_object.decode('utf-8')
            ts = pd.to_datetime(time_object, format='%Y-%m-%d %H-%M-%S')
            np_res_d[i, 0] = (time.mktime(ts.timetuple()))

        IPT_split = eval(config['USER_OPTIONS']['varlist_IPT'])
        for i in range(len(sel_vars)):
            if (sel_vars[i] not in list(map(lambda s: s + '_w', IPT_split))) and \
               (sel_vars[i] not in list(map(lambda s: s + '_T', IPT_split))):
                np_res_d[:, i+1] = results['MET'][sel_vars[i]][...]
            if sel_vars[i] in list(map(lambda s: s + '_w', IPT_split)):
                np_res_d[:, i+1] = results['MET']['qaqc']['IPT_w'][:,IPT_split.index(sel_vars[i][:-2])][...]
            if sel_vars[i] in list(map(lambda s: s + '_T', IPT_split)):
                np_res_d[:, i+1] = results['MET']['qaqc']['IPT_T'][:,IPT_split.index(sel_vars[i][:-2])][...]

        # Some timestamps in the hdf5 files are not in the good format. Might give np_res_d
        # arrays larger than 48.
        # Find the positions of wrong timestamps
        list_del = []
        for timestamp in range(len(np_res_d[:, :])):
            # Not a half-hour
            if np_res_d[timestamp, 0] % (30*60) != 0:
                list_del.append(timestamp)

        # Delete the wrong HH in np_res_d
        for i in reversed(list_del):
            np_res_d = np.delete(np_res_d, i, axis=0)

        return np_res_d


    # =============================================================================
    # %% Read EC-flux, IONICON and GLOVOCS files
    # =============================================================================
    # Store current time
    proc_start_time = datetime.now()
    print('procedure started at ' + proc_start_time.strftime("%d/%m/%Y %H:%M:%S") + '\n')
    
    # Get list of input files
    results_files_list = {'name': [], 'path': [], 'prefix': [], 'date': []}
    cov_files_list = {'name': [], 'path': [], 'prefix': [], 'date': []}
    
    # Load list of results tracer and cov files
    temp_results = list_files_datestring(config['FILES_PATHS']['GEddySoft_results_output_folder'],
                                 config['FILES_NAMES']['GEddySoft_results_files_prefix'],
                                 '.hdf5', 'yyyymmdd',
                                 eval(config['USER_OPTIONS']['date_files_selection']))
    
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        temp_cov = list_files_datestring(config['FILES_PATHS']['GEddySoft_cov_output_folder'],
                                     config['FILES_NAMES']['GEddySoft_cov_files_prefix'],
                                     '.hdf5', 'yyyymmdd',
                                     eval(config['USER_OPTIONS']['date_files_selection']))
    
    for key in results_files_list:
            results_files_list[key] = results_files_list[key] + temp_results[key]
            if config['USER_OPTIONS']['concat_cov'] == 'True':
                cov_files_list[key] = cov_files_list[key] + temp_cov[key]
    
    if results_files_list['name']:
        print('found ' + str(len(results_files_list['name']))
              + ' results files')
        if config['USER_OPTIONS']['concat_cov'] == 'True':
            print('found ' + str(len(cov_files_list['name']))
                  + ' cov files')
            if len(results_files_list['name']) != len(cov_files_list['name']):
                print('Different numbers of results and cov files. Set concat_cov = False'
                      +' or chose other dates.\nProgram interrupted.')
                sys.exit()
    else:
        print('found ' + str(len(results_files_list['name']))
              + ' results files, program interrupted')
        if config['USER_OPTIONS']['concat_cov'] == 'True':
            print('found ' + str(len(cov_files_list['name']))
                  + ' cov files, program interrupted')
        sys.exit()
    
    # Load IONICON product ions table
    os.chdir(os.getcwd()+'\clustering')
    df_IONICON = pd.read_csv(config['FILES_NAMES']['IONICON_file_name'], sep='\s+',
                             names=['mz', 'value', 'chemform', 'comment'])  # list('ABCD')
    df_IONICON['chemform'] = df_IONICON['chemform'] +\
        (' ' + df_IONICON['comment']).fillna('')
    del df_IONICON['comment']
    df_IONICON['chemform'] = df_IONICON['chemform'].map(lambda x: x.lstrip('['))
    df_IONICON['chemform'] = df_IONICON['chemform'].map(lambda x: x.rstrip(']'))
    
    # Load mz_GLOVOCS file
    df_GLOVOCS = pd.read_csv(config['FILES_NAMES']['GLOVOCS_file_name'], sep=',')
    
    
    # =============================================================================
    # %% Concatenate TRACERS, COV and MET
    # =============================================================================
    # To load the data and concatenate tracers, we first use numpy and then switch
    # it to xarrays.
    # For tracers, numpy will be 3D (dim1 = time, dim2 = variables, dim3 = product ions),
    # where measured mz have been attributed to exact mz from IONICON table.
    # For met, numpy will be 2D (dim1 = time, dim2 = variables)
    # For follow_mz, numpy will be 2D (dim1 = time, dim2 = mz labels)
    
    # Additional lines to account for IPT variables to read
    varlist_qaqc = eval(config['USER_OPTIONS']['varlist_qaqc']) 
    # New strings to replace 'IPT'
    IPT_split = eval(config['USER_OPTIONS']['varlist_IPT'])    
    if 'IPT' in varlist_qaqc:
        # split IPT in its different variables
        varlist_qaqc = [new_str for s in varlist_qaqc for new_str in (IPT_split if s == 'IPT' else [s])]   
    varlist_met = eval(config['USER_OPTIONS']['varlist_met']) 
    
    if 'IPT_w' in varlist_met:
        # split IPT in its different variables
        varlist_met = [new_str for s in varlist_met for new_str in
                       (list(map(lambda s: s + '_w', IPT_split)) if s == 'IPT_w' else [s])]
    
    if 'IPT_T' in varlist_met:
        # split IPT in its different variables
        varlist_met = [new_str for s in varlist_met for new_str in
                       (list(map(lambda s: s + '_T', IPT_split)) if s == 'IPT_T' else [s])]
    


    
    # Initialize the tracer 3D multi-day np
    if config['USER_OPTIONS']['concat_option'] == 'IONICON':
        np_res_tracer_md = np.empty((48*len(results_files_list['name']),
                                     len(eval(config['USER_OPTIONS']['varlist_tracer']))
                                     + len(varlist_qaqc) + 2,
                                     len(df_IONICON)))
        if config['USER_OPTIONS']['concat_cov'] == 'True':
            np_res_cov_md = np.empty((48*len(cov_files_list['name']),
                                      int(config['USER_OPTIONS']['n_lags'])+2,
                                      len(df_IONICON)))
    elif config['USER_OPTIONS']['concat_option'] == 'GLOVOCS':
        np_res_tracer_md = np.empty((48*len(results_files_list['name']),
                                     len(eval(config['USER_OPTIONS']['varlist_tracer']))
                                     + len(varlist_qaqc) + 2,
                                     len(df_GLOVOCS)))
        if config['USER_OPTIONS']['concat_cov'] == 'True':
            np_res_cov_md = np.empty((48*len(cov_files_list['name']),
                                      int(config['USER_OPTIONS']['n_lags'])+2,
                                      len(df_GLOVOCS)))
    
    np_res_tracer_md.fill(np.nan)
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        np_res_cov_md.fill(np.nan)
    
    # Initialize the 2D multi-day np
    np_res_met_md = np.empty((48*len(results_files_list['name']),
                              len(varlist_met)+1))
    np_res_met_md.fill(np.nan)
    
    # Initialize the follow_mz 2D multi-day np (used for concatenation of mz)
    follow_mz = np.empty((len(results_files_list['name']), 1400))
    follow_mz.fill(np.nan)
    
    # Fill the 3D multi-day TRACER and COV np with daily outputs, selecting only
    # ion products given by the IONICON/GLOVOCS table
    for d in range(len(results_files_list['name'])):  # loop on daily files
    
        # Load results file into a dictionary
        print('results file ' + results_files_list['name'][d])  # + '\n'
        # Load file only if file_size > 1000 Ko
        file_size_TRACER = os.path.getsize(config['FILES_PATHS']['GEddySoft_results_output_folder']
                                    + '\\' + results_files_list['name'][d])/1000.
        if config['USER_OPTIONS']['results_type'] == 'all':
            # file_size_TRACER_min = 10000
            file_size_TRACER_min = 5000
        elif config['USER_OPTIONS']['results_type'] == 'calibrated':
            file_size_TRACER_min = 500
        if file_size_TRACER > file_size_TRACER_min: 
            results = h5py.File(config['FILES_PATHS']['GEddySoft_results_output_folder']
                                + '\\' + results_files_list['name'][d], mode='r')
        else:
            print('too small')
            continue
        if config['USER_OPTIONS']['concat_cov'] == 'True':    
            file_size_COV = os.path.getsize(config['FILES_PATHS']['GEddySoft_cov_output_folder']
                                    + '\\' + cov_files_list['name'][d])/1000.
            if file_size_COV > 500:
                cov = h5py.File(config['FILES_PATHS']['GEddySoft_cov_output_folder']
                                    + '\\' + cov_files_list['name'][d], mode='r')
            else:
                print('too small')
                continue
    
        # Extract wanted TRACER variables into the 3D np array (time, var, mz)
        results_tracer_d = TRACER_2_3Dnp(results,
                                         eval(config['USER_OPTIONS']['varlist_tracer']),
                                         varlist_qaqc)
        
        # Extract COV functions into the 3D np array (time, lagtime, mz)
        if config['USER_OPTIONS']['concat_cov'] == 'True':
            results_cov_d = COV_2_3Dnp(results, cov, config['USER_OPTIONS']['n_lags'])
    
        # Store all measured mz for this day
        fnonnan = np.argmin(np.where(~np.isnan(results_tracer_d[:, 0, 0]))[0])
        follow_mz[d, 0] = results_tracer_d[fnonnan, 0, 0]
        follow_mz[d, 0:len(results_tracer_d[fnonnan, 0, :])] = results_tracer_d[fnonnan, 0, :]
    
        # Fill np_res_tracer_md (md for "multidays") with:
        # 1) those mz that are within a defined range of mz around the mz
        #  of a known product ion (if concat_option = 'IONICON')
        # 2) those mz that are within a defined range of mz coming from the clustering work
        #  of Bert Verreyken + comparison to GLOVOCS (if concat_option = 'GLOVOCS').
        # Implemented to solve for errors with the deconvolution for different peaks
        # in the peak system)
        if config['USER_OPTIONS']['concat_option'] == 'IONICON':
            # Loop on each mz in the IONICON file
            for j in range(len(df_IONICON)):
                # Compute the mz difference between tracer files and IONICON file
                diff_mz = np.abs(results_tracer_d[0, 1, :] - df_IONICON.loc[j, 'value'])
                # Keep the mz for which the difference is below a treshold
                diff_mz[diff_mz > float(config['USER_OPTIONS']['mz_treshold'])] = np.nan
                if ~np.isnan(diff_mz).all():
                    # Fill-in the multiday file for the mz found in IONICON file
                    idx = np.nanargmin(diff_mz)
                    # Fill-in with results_tracer_d or cov_tracer_d values
                    np_res_tracer_md[d*48:d*48+len(results_tracer_d[:, 0, 0]), :, j]\
                        = results_tracer_d[:, :, idx]
                    if config['USER_OPTIONS']['concat_cov'] == 'True':
                        np_res_cov_md[d*48:d*48+len(results_cov_d[:, 0, 0]), :, j]\
                            = results_cov_d[:, :, idx]
                    # Fill-in with mz value from IONICON file
                    np_res_tracer_md[d*48:d*48+len(results_tracer_d[:, 0, 0]), 1, j]\
                        = df_IONICON.loc[j, 'value']
                    if config['USER_OPTIONS']['concat_cov'] == 'True':
                        np_res_cov_md[d*48:d*48+len(results_cov_d[:, 0, 0]), 1, j]\
                            = df_IONICON.loc[j, 'value']
    
                    # print cases with more than one measured mz in the window around
                    # exact mz, for info purposes at run time
                    nan_nb = np.count_nonzero(~np.isnan(diff_mz))
                    if nan_nb > 1:
                        print(str(nan_nb) +
                              ' mz in the window around searched exact prod. ion for ' +
                              str(df_IONICON.loc[j, 'value']) + ', ')
    
        elif config['USER_OPTIONS']['concat_option'] == 'GLOVOCS':
            # Loop on all the mz for a given day
            for i, mass in enumerate(results_tracer_d[1, 1, :]):
                # Loop on all the mz found during the clustering work
                for j, row in df_GLOVOCS.iterrows():
                    # If a mz from tracer file is found in the GLOVOCS file, write it into
                    # the multiday file
                    if row['mz_min'] <= mass <= row['mz_max']:
                        # Fill-in with results_tracer_d or cov_tracer_d values
                        np_res_tracer_md[d*48:d*48+len(results_tracer_d[:, 0, 0]), :, j]\
                            = results_tracer_d[:, :, i]
                        if config['USER_OPTIONS']['concat_cov'] == 'True':
                            np_res_cov_md[d*48:d*48+len(results_cov_d[:, 0, 0]), :, j]\
                                = results_cov_d[:, :, i]
                        # Fill-in with mz value from GLOVOCS file
                        np_res_tracer_md[d*48:d*48+len(results_tracer_d[:, 0, 0]), 1, j]\
                            = df_GLOVOCS.loc[j, 'mz_avg']
                        if config['USER_OPTIONS']['concat_cov'] == 'True':
                            np_res_cov_md[d*48:d*48+len(results_cov_d[:, 0, 0]), 1, j]\
                                = df_GLOVOCS.loc[j, 'mz_avg']
    
        # Extract wanted MET variables into the 2D np array (time, var)
        results_met_d = MET_2_2Dnp(results, varlist_met)
        
        # Fill np_res_met_md (md for "multidays") with daily values
        np_res_met_md[d*48:d*48+len(results_met_d[:, 0]), :] = results_met_d[:, :]

    # Clean np_res_tracer_md, np_res_cov_md and np_res_met_md by removing empty timestamps
    # (those arrays were created with a default size of 48*ndays but not all days are complete)
    del_list = []
    for n in range(len(np_res_tracer_md[:, 0, 0])):  # loop on time
        if np.isnan(np_res_tracer_md[n, :, :]).all():
            del_list.append(n)
    np_res_tracer_md = np.delete(np_res_tracer_md, del_list, axis=0)
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        np_res_cov_md = np.delete(np_res_cov_md, del_list, axis=0)
    np_res_met_md = np.delete(np_res_met_md, del_list, axis=0)
    
    # Clean np_res_tracer_md, np_res_cov_md and df_IONICON/df_GLOVOCS by removing empty mz
    # (occurence of those products ions never found in results_d)
    del_list = []
    for n in range(len(np_res_tracer_md[0, 0, :])):    # loop on mz
        if np.isnan(np_res_tracer_md[:, :, n]).all():
            del_list.append(n)
    np_res_tracer_md = np.delete(np_res_tracer_md, del_list, axis=2)
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        np_res_cov_md = np.delete(np_res_cov_md, del_list, axis=2)
    if config['USER_OPTIONS']['concat_option'] == 'IONICON':
        df_IONICON.drop(del_list, axis='index', inplace=True)
        df_IONICON.reset_index(inplace=True, drop=True)
    if config['USER_OPTIONS']['concat_option'] == 'GLOVOCS':
        df_GLOVOCS.drop(del_list, axis='index', inplace=True)
        df_GLOVOCS.reset_index(inplace=True, drop=True)
    
    # Fill np_res_tracer_md time and mz cells that are empty because for a part of the time
    # no measured mz fitted the exact mass in IONICON table
    for n in range(len(np_res_tracer_md[:, 0, 0])):  # loop on time
        # index of first non Nan timestamp
        fnonnan = np.argmin(np.where(~np.isnan(np_res_tracer_md[n, 0, :]))[0])
        np_res_tracer_md[n, 0, :] = np_res_tracer_md[n, 0, fnonnan]
        if config['USER_OPTIONS']['concat_cov'] == 'True':
            np_res_cov_md[n, 0, :] = np_res_cov_md[n, 0, fnonnan]
        if config['USER_OPTIONS']['concat_option'] == 'IONICON':
            np_res_tracer_md[n, 1, :] = df_IONICON['value']
            if config['USER_OPTIONS']['concat_cov'] == 'True':
                np_res_cov_md[n, 1, :] = df_IONICON['value']
        if config['USER_OPTIONS']['concat_option'] == 'GLOVOCS':
            np_res_tracer_md[n, 1, :] = df_GLOVOCS['mz_avg']
            if config['USER_OPTIONS']['concat_cov'] == 'True':
                np_res_cov_md[n, 1, :] = df_GLOVOCS['mz_avg']        

    # Remove duplicate values of HH
    unique_indices = np.unique(np_res_tracer_md[:, 0, 0], return_index=True)[1]
    np_res_tracer_md = np_res_tracer_md[unique_indices]
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        np_res_cov_md = np_res_cov_md[unique_indices]
    np_res_met_md = np_res_met_md[unique_indices]
    
    # Find and fill-in with nan the missing HH within missing days
    # Create a list of existing HH
    existing_HH = [datetime.fromtimestamp(ts) for ts in np_res_tracer_md[:, 0, 0]]

    # Find all expected HH
    date_1 = eval(config['USER_OPTIONS']['date_files_selection'])[0][0:17]
    date_2 = eval(config['USER_OPTIONS']['date_files_selection'])[1][0:17]
    date_1 = datetime.strptime(date_1, "%Y_%m_%d__%H_%M")+timedelta(minutes=30)
    date_2 = datetime.strptime(date_2, "%Y_%m_%d__%H_%M")
    if date_2.hour == 0 and date_2.minute == 0:
        date_2 = date_2 + timedelta(days=1)
    current_date = date_1
    expected_HH = []
    while current_date <= date_2:
        expected_HH.append(current_date)
        current_date += timedelta(minutes=30)

    # Transform the lists of expected and existing halfhours into numeric format
    expected_HH = [time.mktime(dt.timetuple()) for dt in expected_HH]
    existing_HH = [time.mktime(dt.timetuple()) for dt in existing_HH]

    # For every missing HH within existing_HH
    for halfhour in expected_HH:
        if halfhour not in existing_HH:
            index = expected_HH.index(halfhour)
        # Create nan lines for tracer, met and cov files
            fill_tracer = np.full((1, np_res_tracer_md.shape[1],
                                   np_res_tracer_md.shape[2]), np.nan)
            fill_tracer[0, 0, :] = halfhour
            for j in range(np_res_tracer_md.shape[2]):
                fill_tracer[:, 1, j] = np_res_tracer_md[:, 1, j][0]
            if config['USER_OPTIONS']['concat_cov'] == 'True':
                fill_cov = np.full((1, np_res_cov_md.shape[1], np_res_cov_md.shape[2]),
                                   np.nan)
                fill_cov[0, 0, :] = halfhour
                for j in range(np_res_cov_md.shape[2]):
                    fill_cov[:, 1, j] = np_res_cov_md[:, 1, j][0]

            fill_met = np.full((1, np_res_met_md.shape[1]), np.nan)
            fill_met[:, 0] = halfhour

            # Add nan lines to tracer, met and cov files
            np_res_tracer_md = np.insert(np_res_tracer_md, index, fill_tracer, axis=0)
            if config['USER_OPTIONS']['concat_cov'] == 'True':
                np_res_cov_md = np.insert(np_res_cov_md, index, fill_cov, axis=0)
            np_res_met_md = np.insert(np_res_met_md, index, fill_met, axis=0)
    
    # Sort the new added HH
    np_res_tracer_md = np_res_tracer_md[np.argsort(np_res_tracer_md[:, 0, 0])]
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        np_res_cov_md = np_res_cov_md[np.argsort(np_res_cov_md[:, 0, 0])]
    np_res_met_md = np_res_met_md[np.argsort(np_res_met_md[:, 0])]

    # # Clean results_md by removing non-requested coumpounds
    # if compoundlist:
    #     del_list = []
    #     for n in range(len(df_IONICON)):    # loop on mz
    #         if not(df_IONICON.iloc[n][2] in compoundlist):
    #             del_list.append(n)
    #     np_res_tracer_md = np.delete(np_res_tracer_md, del_list, axis=2)
    #     df_IONICON.drop(del_list, axis='index', inplace=True)
    #     df_IONICON.reset_index(inplace=True, drop=True)
    
    # =============================================================================
    # %% Manipulate the GLOVOCS IDs for each mz
    # =============================================================================
    if config['USER_OPTIONS']['concat_option'] == 'IONICON':
        mz_chemform = df_IONICON['chemform']
    elif config['USER_OPTIONS']['concat_option'] == 'GLOVOCS':
        # mz as specified in GLOVOCS
        mz_glov = pd.DataFrame()
        # Way of matching clustering and GLOVOCS.
        # ic: in cluster, the m/z value of the matched compound is inside the cluster.
        # nn: nearest neighbor, if no match found within the cluster limits.
        mz_match = pd.DataFrame()
        # Chemical formula
        mz_chemform = pd.DataFrame()
        # Compound name
        mz_name = pd.DataFrame()
    
        # In df_GLOVOCS, mz information is separated by '_'. Create a dataframe for each kind
        # of information.
        for column in df_GLOVOCS.loc[:, 'ID0':'ID24']:
            if df_GLOVOCS.loc[:, column].isnull().values.all():
                continue
            else:
                mz_glov[column] = df_GLOVOCS.loc[:, column].str.split('_', expand=True)\
                    .loc[:, 0]
                mz_match[column] = df_GLOVOCS.loc[:, column].str.split('_', expand=True)\
                    .loc[:, 1]
                mz_chemform[column] = df_GLOVOCS.loc[:, column].str.split('_', expand=True)\
                    .loc[:, 2]
                mz_name[column] = df_GLOVOCS.loc[:, column].str.split('_', expand=True)\
                    .loc[:, 3]
    
        # In the xarray Dataset, we cannot put in a dataframe. For each kind of mz
        # information, the columns are joined and separated by '/'.
        mz_glov = np.array(['/'.join(row) for row in
                            mz_glov.loc[:, :].to_numpy().astype(str)])
        mz_match = np.array(['/'.join(row) for row in
                             mz_match.loc[:, :].to_numpy().astype(str)])
        mz_chemform = np.array(['/'.join(row) for row in
                                mz_chemform.loc[:, :].to_numpy().astype(str)])
        mz_name = np.array(['/'.join(row) for row in
                            mz_name.loc[:, :].to_numpy().astype(str)])
    
    # =============================================================================
    # %% Read the 3D np into an xarray
    # =============================================================================
    # Create the 2 dimensions of the xarray (time, mz)
    dt = [datetime.fromtimestamp(s) for s in np_res_tracer_md[:, 0, 0]]
    mz = np_res_tracer_md[0, 1, :]
    # chemform = df_IONICON.chemform
    
    # Create the 3d dimension of the cov xarray
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        lag = np.linspace(0, int(config['USER_OPTIONS']['n_lags'])-1,
                          int(config['USER_OPTIONS']['n_lags']))
    
    # Remove time and mz columns to only keep measured values
    np_res_tracer_md_values = np_res_tracer_md[:, 2:, :]
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        np_res_cov_md_values = np_res_cov_md[:, 2:, :]
    
    # Load the results data in a dict for each tracer variable, from the 3D np
    dict_xr_results = {}
    
    for i in range(len(eval(config['USER_OPTIONS']['varlist_tracer']))):
        # the 2 dimensions must be stated in the dict
        dict_xr_results[eval(config['USER_OPTIONS']['varlist_tracer'])[i]] =\
            xr.Variable(('time', 'mz'),
                        np_res_tracer_md_values[:, i, :].reshape
                        (np_res_tracer_md_values.shape[0],
                         np_res_tracer_md_values.shape[2]),
                        # Attributes of each variable are also added
                        {'name': results['IRGA']['0'][eval(config['USER_OPTIONS']
                                                            ['varlist_tracer'])[i]].attrs['name'],
                         'units': results['IRGA']['0'][eval(config['USER_OPTIONS']
                                                              ['varlist_tracer'])[i]].attrs['units']})
    
    for i in range(len(varlist_qaqc)):
        # the 2 dimensions must be stated in the dict
        if varlist_qaqc[i] not in IPT_split:
            dict_xr_results[varlist_qaqc[i]] =\
                xr.Variable(('time', 'mz'),
                            np_res_tracer_md_values[:, i+len(eval(config['USER_OPTIONS']
                                                                 ['varlist_tracer'])), :]
                            .reshape(np_res_tracer_md_values.shape[0],
                                     np_res_tracer_md_values.shape[2]),
                            # Attributes of each variable are also added
                            {'name': results['IRGA']['0']['qaqc'][varlist_qaqc[i]].attrs['name'],
                             'units': results['IRGA']['0']['qaqc'][varlist_qaqc[i]].attrs['units']})
        else:
            # for IPT variables no attributes are available
            dict_xr_results[varlist_qaqc[i]] =\
                xr.Variable(('time', 'mz'),
                            np_res_tracer_md_values[:, i+len(eval(config['USER_OPTIONS']
                                                                 ['varlist_tracer'])), :]
                            .reshape(np_res_tracer_md_values.shape[0],
                                     np_res_tracer_md_values.shape[2]),
                            # Attributes of each variable are also added
                            {'name': varlist_qaqc[i],
                             'units': '%'})

    # Add in the dict the met data
    for i in range(len(varlist_met)):
        # the 2 dimensions must be stated in the dict
        if (varlist_met[i] not in list(map(lambda s: s + '_w', IPT_split))) and (varlist_met[i] not in list(map(lambda s: s + '_T', IPT_split))):
            if 'name' in results['MET'][varlist_met[i]].attrs:
                dict_xr_results[varlist_met[i]] =\
                    xr.Variable('time',
                                np_res_met_md[:, 1+i],
                                # Attributes of each variable are also added
                                {'name': results['MET'][varlist_met[i]].attrs['name'],
                                 'units': results['MET'][varlist_met[i]].attrs['units']})
            else:
                print('Name and units missing for var ' + varlist_met[i])
                dict_xr_results[varlist_met[i]] =\
                    xr.Variable('time',
                                np_res_met_md[:, 1+i],
                                # Attributes of each variable are also added
                                {'name': 'name missing',
                                 'units': 'units missing'})                
        else:
            dict_xr_results[varlist_met[i]] =\
                xr.Variable('time',
                            np_res_met_md[:, 1+i],
                            # Attributes of each variable are also added
                            {'name': varlist_met[i],
                             'units': '%'})
  
    
    # Add in the results dict the available information regarding mz
    if config['USER_OPTIONS']['concat_option'] == 'IONICON':
        dict_xr_results['mz_chemform'] = xr.Variable('mz', mz_chemform)
    elif config['USER_OPTIONS']['concat_option'] == 'GLOVOCS':
        dict_xr_results['mz_glov'] = xr.Variable('mz', mz_glov)
        dict_xr_results['mz_match'] = xr.Variable('mz', mz_match)
        dict_xr_results['mz_chemform'] = xr.Variable('mz', mz_chemform)
        dict_xr_results['mz_name'] = xr.Variable('mz', mz_name)
    
    # Create the xarray Dataset from the dict and assign the 2 coordinates
    ds_results = (xr.Dataset(dict_xr_results,
                     coords={'time': dt,
                             'mz': mz}))

    # Output csv if requested
    if config['USER_OPTIONS']['output_csv'] == 'True':
        csv_folder = config['FILES_PATHS']['concat_results_output_folder']\
            + '\\csv_'+ eval(config['USER_OPTIONS']['date_files_selection'])[0][0:10]\
                + '_'+eval(config['USER_OPTIONS']['date_files_selection'])[1][5:10]\
                    + '_'+ config['USER_OPTIONS']['concat_option']
        # Verify if the csv folder exists
        if not os.path.exists(csv_folder):
            os.makedirs(csv_folder)
        os.chdir(csv_folder)
        data = {}
        for var in varlist_met:
            data[var + ' ' + str(ds_results[var].attrs)]=ds_results[var].values
            
        df = pd.DataFrame(data, index=ds_results['time'].values)
        df.index.name = 'TIMESTAMP_END'
        df.to_csv(config['USER_OPTIONS']['meas_level'] + '_MET.csv', index=True)

        col_title = ['mz_' + element for element in ds_results['mz'].values.astype(str)]
        for var in eval(config['USER_OPTIONS']['varlist_tracer']):
            col_title_var = [element + '_' + str(ds_results[var].attrs)
                             for element in col_title]
            df = pd.DataFrame(ds_results[var].values, index=ds_results['time'].values,
                              columns=col_title_var)
            df.index.name = 'TIMESTAMP_END'
            df.to_csv(config['USER_OPTIONS']['meas_level'] + '_TRACER_' + var +'.csv',
                      index=True)
        for var in varlist_qaqc:
            col_title_var = [element + '_' + str(ds_results[var].attrs)
                             for element in col_title]
            df = pd.DataFrame(ds_results[var].values, index=ds_results['time'].values,
                              columns=col_title_var)
            df.index.name = 'TIMESTAMP_END'
            df.to_csv(config['USER_OPTIONS']['meas_level'] + '_QAQC_' + var +'.csv',
                      index=True)    
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        # Load the cov data in a dict from the 3D np
        dict_xr_cov = {}
        dict_xr_cov['cov_func'] = xr.Variable(('time', 'lag', 'mz'),
                                            np_res_cov_md_values)
        
        # Add in the cov dict the available information regarding mz
        if config['USER_OPTIONS']['concat_option'] == 'IONICON':
            dict_xr_cov['mz_chemform'] = xr.Variable('mz', mz_chemform)
        elif config['USER_OPTIONS']['concat_option'] == 'GLOVOCS':
            dict_xr_cov['mz_glov'] = xr.Variable('mz', mz_glov)
            dict_xr_cov['mz_match'] = xr.Variable('mz', mz_match)
            dict_xr_cov['mz_chemform'] = xr.Variable('mz', mz_chemform)
            dict_xr_cov['mz_name'] = xr.Variable('mz', mz_name)
        
        # Create the xarray cov Dataset from the dict and assign the 3 coordinates
        ds_cov = (xr.Dataset(dict_xr_cov,
                         coords={'time': dt,
                                 'mz': mz,
                                 'lag': lag}))
    
    # Save the xarrays in netCDF format
    os.chdir(config['FILES_PATHS']['concat_results_output_folder'])
    ds_results.to_netcdf('results_' + eval(config['USER_OPTIONS']['date_files_selection'])[0][0:10]
                  + '_'+eval(config['USER_OPTIONS']['date_files_selection'])[1][5:10] + '_'
                  + config['USER_OPTIONS']['concat_option']
                  + '.nc')
    
    if config['USER_OPTIONS']['concat_cov'] == 'True':
        os.chdir(config['FILES_PATHS']['concat_cov_output_folder'])
        ds_cov.to_netcdf('cov_' + eval(config['USER_OPTIONS']['date_files_selection'])[0][0:10]
                      + '_'+eval(config['USER_OPTIONS']['date_files_selection'])[1][5:10] + '_'
                      + config['USER_OPTIONS']['concat_option']
                      + '.nc')
