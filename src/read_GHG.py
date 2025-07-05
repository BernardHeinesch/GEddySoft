"""Functions for reading and processing LI-COR GHG high-frequency data files.

This module provides functionality to read and process high-frequency data from LI-COR GHG files,
specifically focusing on the SMARTFLUX system output. It includes tools for:

* Reading and extracting data from zipped GHG files
* Processing diagnostic values from LI-7200 gas analyzer
* Handling AGC (Automatic Gain Control) values and other diagnostic flags

Author: Ariane Faures
Created: October 5, 2021
"""

import zipfile
import glob
import pandas as pd
import os
import numpy as np

# %% List of functions
def read_GHG (raw_file, raw_format='ghg', unzip_path=None):
    """Read and extract high-frequency data from LI-COR SMARTFLUX GHG files.

    This function handles the reading of high-frequency eddy covariance data from
    LI-COR SMARTFLUX GHG files. It extracts both data and metadata from zipped GHG
    files and returns them as pandas DataFrames.

    Parameters
    ----------
    raw_file : str
        Path to the GHG file to process
    raw_format : str, optional
        Format of the raw data file, currently only 'ghg' is supported
    unzip_path : str, optional
        Directory where the GHG file should be temporarily extracted. If None,
        uses the same directory as the GHG file

    Returns
    -------
    list
        A list containing:

        - file_header : pandas.DataFrame
          Header information from the data file (first 6 lines)
        - file_data : pandas.DataFrame
          High frequency data with variable names as columns
        - data_name : str
          Path to the extracted data file
        - metadata_name : str
          Path to the extracted metadata file

    Notes
    -----
    The function automatically cleans up extracted files after reading them.
    """

    # Read the header (first 7 lines of the file)
    # And afterwards read the many body (with variable names as columns)
    if raw_format == 'ghg':
        with zipfile.ZipFile(raw_file, 'r') as zip_ref:
            zip_ref.extractall(unzip_path)

        data_name = glob.glob(unzip_path + '/' +'*.data')
        metadata_name = glob.glob(unzip_path + '/' + '*.metadata')

        with open(data_name[-1], mode='r') as file:
            file_header = pd.read_table(file, nrows = 6, header = None)
        with open(data_name[-1], mode='r') as file:    
            file_data = pd.read_table(file, header = 7)

    file.close()

    os.remove(data_name[0])
    os.remove(metadata_name[0])    

    return([file_header,file_data,data_name, metadata_name])



# Legacy version of read_diag_val with additional debugging output
# Currently unused - the production version is defined below
# def read_diag_val(data, data_name_short):
#     """Process diagnostic values from LI-7200 gas analyzer data.
# 
#     This function processes the diagnostic values from LI-7200 gas analyzer data,
#     converting binary diagnostic flags into meaningful status indicators for various
#     instrument components.
# 
#     Parameters
#     ----------
#     data : pandas.DataFrame
#         DataFrame containing the raw data with a 'Diagnostic Value' column
#     data_name_short : str
#         Short identifier for the data file being processed
# 
#     Returns
#     -------
#     pandas.DataFrame
#         DataFrame containing diagnostic counts/values for:
#         - AGC (Automatic Gain Control)
#         - Sync
#         - PLL
#         - Detector
#         - Chopper
#         - DeltaPressure
#         - Aux_input
#         - Tinlet
#         - Toutlet
#         - Head detect
#         - Anemometer Diagnostics
# 
#     Notes
#     -----
#     AGC values are compared against a reference of 100.05 and stored as means.
#     Other diagnostic flags are counted when they indicate an issue (value != 1).
#     """
#     print('Reading diagnostic data of files')
#     diag_val = data.loc[:,'Diagnostic Value'].copy() # Copy of the initial column: change!!!
# 
#     print('Starting to count the number of times a flag was raised')
# 
#     def int_to_binary_vectorized(arr):
#         """Convert integer array to binary strings with leading zeros.
# 
#         Parameters
#         ----------
#         arr : numpy.ndarray
#             Array of integers to convert
# 
#         Returns
#         -------
#         numpy.ndarray
#             Array of binary strings, each padded to 7 digits
#         """
#         int2binary = np.vectorize(lambda x: '000' + format(x, 'b'))
#         return int2binary(arr)
#     
#     # Example usage
#     arr = diag_val.to_numpy()
#     diag_val = int_to_binary_vectorized(arr)
#     
#     
#     def split_bin_diag_vectorized(diag_val,data_name_short):
#         """Process binary diagnostic values into component-specific flags.
# 
#         Parameters
#         ----------
#         diag_val : numpy.ndarray
#             Array of binary strings representing diagnostic values
#         data_name_short : str
#             Short identifier for the data file being processed
# 
#         Returns
#         -------
#         pandas.DataFrame
#             DataFrame containing diagnostic counts for each component
#         """
#         # Create the binary array from the input series
#         binary_array = np.array([list(x) for x in diag_val.astype(str)])
#         
#         # Create a mask for the AGC columns: 4 last bits, so columns of the array
#         agc_mask = binary_array[:,-4:].astype(int)
#         num_rows, num_cols = agc_mask.shape
#         # Use the bitwise left shit operator (<<) to comnbine the binary digits
#         # and get one number per row without looping over all the records
#         result = np.zeros(num_rows, dtype=np.int)
#         for i in range(num_cols):
#             result = result << 1
#             result = result | agc_mask[:, i]
#             
#         # Convert the agc diagnostic into its real value to compare it to the
#         # reference value (*6.67 according to the manual)
#         agc_column = (result * 6.67).astype(float) 
#         # Count the number of times it isn't the reference value (100.05) with
#         # a limit of tolerance
#         # agc_count = (~np.isclose(agc_column, 100.05)).sum()
#         agc_count = np.mean(agc_column) # put the mean value instead of the count
#         
#         # Create a mask for the remaining columns
#         other_mask = binary_array[:,3:-4].astype(int)
#         other_count = (~(other_mask == 1)).sum(axis=0)
#         
#         # Create the count dataframe from the counts. First are the agc counts
#         # and then the other counts but flipped since the bits must be read
#         # right to left to correspond to the variable names
#         diag_count = pd.DataFrame([agc_count] + list(np.flip(other_count)),
#                                   columns = [data_name_short[:-9]])
#         diag_count = diag_count.T
#         diag_count.columns = ['AGC', 'Sync', 'PLL', 'Detector', 
#                                            'Chopper', 'DeltaPressure',
#                                            'Aux_input', 'Tinlet', 'Toutlet', 
#                                            'Head detect']
#         
#         return diag_count
# 
#     diag_count = split_bin_diag_vectorized(diag_val,data_name_short)
#     diag_count['Anemometer Diagnostics'] = data.loc[:,'Anemometer Diagnostics'].sum()
#         
#         # print_progress(k, diag_val)
#     agc = pd.Series(((diag_count.iloc[:,0])<100).any())
#     agc.index= ['AGC']
#     diag_flag = pd.concat([agc,((diag_count.iloc[:,1:])!=0).any()])
#     if not((diag_count.loc[:,diag_flag]).empty):     
#         print('File ' + data_name_short)
#         print('Number and type of LI7200 diagnostic flags raised:')
#         print(diag_count.loc[:,diag_flag])
#    
#     return(diag_count)







def read_diag_val(data):
    """Process LI-7200 gas analyzer diagnostic values.

    This function processes diagnostic values from a LI-7200 gas analyzer,
    converting binary diagnostic flags into meaningful status indicators.
    Each diagnostic value is a binary number where specific bits indicate
    the status of different analyzer components.

    Parameters
    ----------
    data : pandas.Series or numpy.ndarray
        Series or array of diagnostic values from the analyzer

    Returns
    -------
    pandas.DataFrame
        DataFrame containing diagnostic information for each component:

        AGC : float
            Mean AGC value (Automatic Gain Control, scaled by 6.67)
        Sync : int
            Count of synchronization issues
        PLL : int
            Count of Phase-Locked Loop issues
        Detector : int
            Count of detector issues
        Chopper : int
            Count of chopper wheel issues
        DeltaPressure : int
            Count of pressure differential issues
        Aux_input : int
            Count of auxiliary input issues
        Tinlet : int
            Count of inlet temperature issues
        Toutlet : int
            Count of outlet temperature issues
        Head detect : int
            Count of head detection issues
    """

    diag_val = data.astype(int)

    def int_to_binary_vectorized(arr):
        """Convert integer array to binary strings with leading zeros.

        Parameters
        ----------
        arr : numpy.ndarray
            Array of integers to convert

        Returns
        -------
        numpy.ndarray
            Array of binary strings, each padded with leading zeros
        """
        int2binary = np.vectorize(lambda x: '000' + format(x, 'b'))
        return int2binary(arr)

    diag_val = int_to_binary_vectorized(diag_val)

    def split_bin_diag_vectorized(diag_val):
        """Process binary diagnostic values into component-specific flags.

        Parameters
        ----------
        diag_val : numpy.ndarray
            Array of binary strings representing diagnostic values

        Returns
        -------
        pandas.DataFrame
            DataFrame containing diagnostic counts/means for each component
        """
        # Create the binary array from the input series
        binary_array = np.array([list(x) for x in diag_val.astype(str)])
        
        # Create a mask for the AGC columns: 4 last bits, so columns of the array
        agc_mask = binary_array[:,-4:].astype(int)
        num_rows, num_cols = agc_mask.shape
        # Use the bitwise left shit operator (<<) to comnbine the binary digits
        # and get one number per row without looping over all the records
        result = np.zeros(num_rows, dtype=int)
        for i in range(num_cols):
            result = result << 1
            result = result | agc_mask[:, i]
            
        # Convert the agc diagnostic into its real value to compare it to the
        # reference value (*6.67 according to the manual)
        agc_column = (result * 6.67).astype(float) 
        # store the mean value
        agc_count = np.mean(agc_column)
        
        # Create a mask for the remaining columns
        other_mask = binary_array[:,3:-4].astype(int)
        other_count = (~(other_mask == 1)).sum(axis=0)
        
        # Create the count dataframe from the counts. First are the agc counts
        # and then the other counts but flipped since the bits must be read
        # right to left to correspond to the variable names
        diag_count = pd.DataFrame([agc_count] + list(np.flip(other_count)))  #,columns = [data_name_short[:-9]]
        diag_count = diag_count.T
        diag_count.columns = ['AGC', 'Sync', 'PLL', 'Detector', 
                                           'Chopper', 'DeltaPressure',
                                           'Aux_input', 'Tinlet', 'Toutlet', 
                                           'Head detect']
        
        return diag_count

    diag_count = split_bin_diag_vectorized(diag_val)

    # print('Number and type of LI7200 diagnostic flags raised:')
    # print(diag_count.loc[:,(diag_count!=0).any()])
   
    return(diag_count)
    

