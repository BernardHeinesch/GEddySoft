"""
Starting routine to run GEddySoft.
Reads the ini file and allows to run the code either in non-multiprocess mode
or in multi-process mode (one process per day of input data)
"""

import multiprocessing as mp
import datetime
import os
import sys
from GEddySoft_main import GEddySoft_main
from get_ini import get_ini
from get_list_input_files import get_list_input_files
from check_missing_dates_parallel_processing import check_missing_dates_parallel_processing
from map_sonic2tracer import map_sonic2tracer


# %% parallel processing routine
def run_pool(days_to_process, log_filename):
    with mp.Pool(cpu) as pool:
        try:
            # Save original stdout before redirecting
            original_stdout = sys.stdout
            # Redirect stdout to devnull for all processes
            sys.stdout = open(os.devnull, 'w')
            results = pool.starmap(GEddySoft_main, [(day, ini, log_filename) for day in days_to_process])
            # Restore stdout
            sys.stdout = sys.__stdout__
            print(f'processed days: {results}' + '\n')
            OF.write(f'successfully processed days: {results}\n')
        except Exception as e:
            # Make sure to restore stdout even if there's an error
            sys.stdout = original_stdout
            print(f'error during processing: {str(e)}')
            OF.write(f'error during processing: {str(e)}\n')
        finally:
            pool.close()
            pool.join()
            print('pool closed')


# %% main
if __name__ == '__main__':

    # -------------------------------------------------------------------------------------------------------
    # user choices
    ini_filename = r'../examples/metadata/GEddySoft_parameters_VOC_top.ini'

    # or read it from command line
    if len(sys.argv) > 1:
        ini_filename = sys.argv[1]
    # -------------------------------------------------------------------------------------------------------

    # store current time
    proc_start_time = datetime.datetime.now()
    print('run started at ' + proc_start_time.strftime("%d/%m/%Y %H:%M:%S") + '\n')

    # get ini information
    ini = get_ini(ini_filename)

    # Open logfile in write mode, and close it afterwards
    # start_timestamp = proc_start_time.strftime("%Y%m%d_%H%M%S")
    # log_filename = f"{ini['files']['log_filepath']}/logfile_{start_timestamp}.csv"
    log_filename = f"{ini['files']['log_filepath']}/logfile.csv"
    OF = open(log_filename, 'w')
    OF.write('\n--- New Processing Session ---\n\n')
    OF.write('run started at ' + proc_start_time.strftime("%d/%m/%Y %H:%M:%S") + '\n' + '\n')
    print('*************** reading metadata ***************' + '\n')
    OF.write('*************** reading metadata ***************' + '\n' + '\n')
    OF.flush()
    OF.close()

    # Open logfile in append mode instead of write mode
    OF = open(log_filename, 'a')

    if ini['run_param']['MULTIPROCESSING']:

        # run in multiprocessing mode

        # get list of sonic and tracer input files
        all_sonic_files_list, all_tracer_files_list = get_list_input_files(ini)

        # map sonic file list to tracer file list if asked for
        if ini['run_param']['MAP_SONIC2TRACER'] and ini['run_param']['CONCENTRATION_TYPE']:
            all_sonic_files_list = map_sonic2tracer(all_sonic_files_list, all_tracer_files_list)
            print('sonic files mapped to tracer files\n')

        print('*************** processing data ***************' + '\n')

        # get uniques yyyy_mm_dd in the sonic file list
        unique_days = list(set(list(map(lambda x: x[len(ini['files']['sonic_files_prefix']):len(ini['files']['sonic_files_prefix']) + 10], all_sonic_files_list['name']))))
        unique_days.sort()

        unique_days = [date.replace('-', '_') for date in unique_days]

        cpu = mp.cpu_count()
        print('running in multithread mode with ' + str(cpu) + ' CPUs... (console muted)\n')
        OF.write(f'running in multithread mode with {cpu} CPUs\n')

        # First attempt with all days
        run_pool(unique_days, log_filename)

        # Check for missing dates
        missing_dates = check_missing_dates_parallel_processing(
            ini['files']['output_folder'],
            ini['files']['output_files_prefix'],
            unique_days
        )

        if missing_dates:
            print(f"First pass: Found {len(missing_dates)} missing dates:")
            OF.write(f"First pass: Found {len(missing_dates)} missing dates:\n")
            print('[')
            for date in sorted(missing_dates): print(f"{date},")
            print(']')

            print("Retrying missing dates...")
            OF.write("Retrying missing dates...\n")
            # Second attempt with missing dates only
            run_pool(missing_dates, log_filename)

            # Check again after retry
            final_missing_dates = check_missing_dates_parallel_processing(
                ini['files']['output_folder'],
                ini['files']['output_files_prefix'],
                unique_days
            )

            if final_missing_dates:
                print(f"WARNING: After retry, still found {len(final_missing_dates)} unprocessed dates:")
                OF.write(f"WARNING: After retry, still found {len(final_missing_dates)} unprocessed dates:\n")
                for date in sorted(final_missing_dates): print(f"  - {date}")
                print("Consider investigating these dates manually.")
            else:
                print("All dates successfully processed after retry!")
                OF.write("All dates successfully processed after retry!\n")
        else:
            print("No missing dates found - all files were processed successfully!")

    else:

        # run in normal mode
        GEddySoft_main('no_multithread', ini, log_filename)
        
    # display time needed
    proc_end_time = datetime.datetime.now()
    print('finished at ' + proc_end_time.strftime("%d/%m/%Y %H:%M:%S"))
    print('run took ' + str((proc_end_time - proc_start_time)))
    OF.write('finished at ' + proc_end_time.strftime("%d/%m/%Y %H:%M:%S") + "\n")
    OF.write('. run took ' + str((proc_end_time - proc_start_time)) + "\n")

    # permanent closure of the log file
    OF.close()
