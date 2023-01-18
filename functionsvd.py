import numpy as np
from tqdm import tqdm
import pickle
from time import perf_counter
from scipy.sparse.linalg import svds
import datetime
import h5py
import glob
from scipy.signal import detrend
from numpy.fft import fftshift, fft2, fftfreq
from datetime import datetime as DT
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def data_wrangler(cable,record_length,t0):
    if cable == 'seadasn':
        prefix = 'seadasn'
        network_name = 'SeaDAS-N'
        if t0 < datetime.datetime(2022, 6, 20, 0, 0, 0):
            datastore='/data/data0/seadasn_2022-03-17_2022-06-20/'
        elif (t0 >= datetime.datetime(2022, 6, 20, 0, 0, 0)) and (t0 < datetime.datetime(2022, 10, 7, 0, 0, 0)):
            datastore='/data/data7/seadasn_2022-06-21_2022-10-06/'
        else:
            datastore='/data/data3/seadasn/'

    elif cable == 'whidbey':
        prefix = 'whidbey'
        network_name='Whidbey-DAS'
        if t0 < datetime.datetime(2022,10,23,4,50,0):
            datastore = '/data/data5/Converted/'
        else:
            datastore = '/data/data6/whidbey/'
        
    return prefix, network_name, datastore

def sintela_to_datetime(sintela_times):
    '''
    returns an array of datetime.datetime 
    ''' 
    days1970 = datetime.date(1970, 1, 1).toordinal()

    # Vectorize everything
    converttime = np.vectorize(datetime.datetime.fromordinal)
    addday_lambda = lambda x : datetime.timedelta(days=x)
    adddays = np.vectorize(addday_lambda )
    
    day = days1970 + sintela_times/1e6/60/60/24
    thisDateTime = converttime(np.floor(day).astype(int))
    dayFraction = day-np.floor(day)
    thisDateTime = thisDateTime + adddays(dayFraction)

    return thisDateTime



def open_sintela_file(file_base_name,t0,pth,
                      chan_min=0,
                      chan_max=-1,
                      number_of_files=1,
                      verbose=False,
                      pad=False):

    data = np.array([])
    time = np.array([])

    
    dt = datetime.timedelta(minutes=1) # Assume one minute file duration
    this_files_date = t0
    
    for i in range(number_of_files):
        
        # Construct the "date string" part of the filename
        date_str = this_files_date.strftime("%Y-%m-%d_%H-%M")
    
        # Construct the PARTIAL file name (path and name, but no second or filenumber):
#         this_file = f'{pth}{file_base_name}_{date_str}_UTC_{file_number:06}.h5'
        partial_file_name = f'{pth}{file_base_name}_{date_str}'
        file_search = glob.glob(f'{partial_file_name}*h5')
        if verbose:
            print(f'Searching for files matching: {partial_file_name}*h5')
        if len(file_search) > 1:
            raise ValueError("Why are there more than one files? That shouldn't be possible!")
        elif len(file_search) == 0:
            raise ValueError("Why are there ZERO files? That shouldn't be possible!")
        else:
            this_file = file_search[0]
        
        try:
            f = h5py.File(this_file,'r')
            this_data = np.array(
                f['Acquisition/Raw[0]/RawData'][:,chan_min:chan_max])
            this_time = np.array(
                f['Acquisition/Raw[0]/RawDataTime'])
            
            if i == 0:
                time = sintela_to_datetime(this_time)
                data = this_data
                attrs=dict(f['Acquisition'].attrs)
            else:
                data = np.concatenate((data, this_data ))
                time = np.concatenate((time, this_time ))
                
        except Exception as e: 
            print('File problem with: %s'%this_file)
            print(e)
           


        this_files_date = this_files_date + dt
    
    #if pad==True:
        # Add columns of zeros to give data matrix the correct dimensions
        
    return data, time, attrs