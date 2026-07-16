import zipfile
import glob
import pandas as pd
import os
import numpy as np

# %% List of functions
def read_GHG (raw_file, raw_format='ghg', unzip_path=None):
    """
    Functions to read and visualise HF data from 
        - .GHG (LICOR SMARTFLUX)
    
    Generalisation to all Raw files. unzip path set to None as default because 
    needed only for GHG files
    List files in the input directory and unzip them in the unzip directory
    Read the data and metadata files extracted   
    Inputs:
        - specific file name
        - directory of the unzipped files
    Outputs:
        - pandas dataframe with high frequency data
        - pandas dataframe with header of data file
        - pandas dataframe with metadata information
        - unzipped data file name
        - unzipped metadata file name

    comments
    --------
    Created on Tue Oct  5 15:08:00 2021
    @author: Ariane FAURES
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



def read_diag_val(data, data_name_short):
    print('Reading diagnostic data of files')
    diag_val = data.loc[:,'Diagnostic Value'].copy() # Copy of the initial column: change!!!

    print('Starting to count the number of times a flag was raised')

    def int_to_binary_vectorized(arr):
        int2binary = np.vectorize(lambda x: '000' + format(x, 'b'))
        return int2binary(arr)
    
    # Example usage
    arr = diag_val.to_numpy()
    diag_val = int_to_binary_vectorized(arr)
    
    
    def split_bin_diag_vectorized(diag_val,data_name_short):
    # Create the binary array from the input series
        binary_array = np.array([list(x) for x in diag_val.astype(str)])
        
        # Create a mask for the AGC columns: 4 last bits, so columns of the array
        agc_mask = binary_array[:,-4:].astype(int)
        num_rows, num_cols = agc_mask.shape
        # Use the bitwise left shit operator (<<) to comnbine the binary digits
        # and get one number per row without looping over all the records
        result = np.zeros(num_rows, dtype=np.int)
        for i in range(num_cols):
            result = result << 1
            result = result | agc_mask[:, i]
            
        # Convert the agc diagnostic into its real value to compare it to the
        # reference value (*6.67 according to the manual)
        agc_column = (result * 6.67).astype(float) 
        # Count the number of times it isn't the reference value (100.05) with
        # a limit of tolerance
        # agc_count = (~np.isclose(agc_column, 100.05)).sum()
        agc_count = np.mean(agc_column) # put the mean value instead of the count
        
        # Create a mask for the remaining columns
        other_mask = binary_array[:,3:-4].astype(int)
        other_count = (~(other_mask == 1)).sum(axis=0)
        
        # Create the count dataframe from the counts. First are the agc counts
        # and then the other counts but flipped since the bits must be read
        # right to left to correspond to the variable names
        diag_count = pd.DataFrame([agc_count] + list(np.flip(other_count)),
                                  columns = [data_name_short[:-9]])
        diag_count = diag_count.T
        diag_count.columns = ['AGC', 'Sync', 'PLL', 'Detector', 
                                           'Chopper', 'DeltaPressure',
                                           'Aux_input', 'Tinlet', 'Toutlet', 
                                           'Head detect']
        
        return diag_count

    diag_count = split_bin_diag_vectorized(diag_val,data_name_short)
    diag_count['Anemometer Diagnostics'] = data.loc[:,'Anemometer Diagnostics'].sum()
        
        # print_progress(k, diag_val)
    agc = pd.Series(((diag_count.iloc[:,0])<100).any())
    agc.index= ['AGC']
    diag_flag = pd.concat([agc,((diag_count.iloc[:,1:])!=0).any()])
    if not((diag_count.loc[:,diag_flag]).empty):     
        print('File ' + data_name_short)
        print('Number and type of LI7200 diagnostic flags raised:')
        print(diag_count.loc[:,diag_flag])
   
    return(diag_count)







def read_diag_val(data):
    """
    Read diagnostic value and flag the malfunctioning instrument part
    Inputs:
        - pandas dataframe with variable names in column name
    Outputs:
        - Diagnostic value translated into flag for each instrument part
    """

    diag_val = data.astype(int)

    def int_to_binary_vectorized(arr):
        int2binary = np.vectorize(lambda x: '000' + format(x, 'b'))
        return int2binary(arr)

    diag_val = int_to_binary_vectorized(diag_val)

    def split_bin_diag_vectorized(diag_val):
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
    

