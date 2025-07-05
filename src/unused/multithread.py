from GEddySoft_TOF_step1 import GEddySoft
from multiprocess import Pool
import datetime

# %%


def running(i):
    if len(str(i)) == 1:
        print('0'+str(i))
        GEddySoft('0'+str(i))
    else:
        print(str(i))
        GEddySoft(str(i))
    return None

# user choices
multiprocessing = True

# define the days to be processed ([0] if the info in the ini should be used instead)
# must be an iterable (e.g. list or range)
days = [6,11]  # range(1,32)

if multiprocessing:

    if __name__ == '__main__':
    
        # OutputFile= args.Output_File
        OF = open('D:\Documents\Science\Programmation\Projets Python\GEddySoft_TOF\GEddySoft_TOF_v3.0\output\logfile.csv', 'w')
    
        # store current time
        proc_start_time = datetime.datetime.now()
        OF.write('procedure started at ' + proc_start_time.strftime("%d/%m/%Y %H:%M:%S") + "\n")
    
        # # set processes number to be used
        # pool = Pool(7)
    
        # print('running...')
    
        # # call multiprocessing
        # ResultList = pool.map(running, days)
        # # trial 2
        # #ResultList = pool.starmap(running, range(1,2))
        # #pool.close()
        # #pool.join()
        
        import multiprocessing
        num_processes = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(num_processes)
        results = pool.map_async(running, days)
        pool.close()
        pool.join()
    
    
        # display time needed
        proc_end_time = datetime.datetime.now()
        OF.write('finished at ' + proc_end_time.strftime("%d/%m/%Y %H:%M:%S") + '. Procedure took ' + str((proc_end_time - proc_start_time)) + "\n")
        OF.close
    
        print('finished')

else:

    GEddySoft('00')



