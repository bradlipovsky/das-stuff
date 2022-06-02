import numpy as np
import datetime
import h5py

import glob
def get_file_number(pth,prefix,t0):
    
    datestr = '{d.year}-{d.month:02}-{d.day:02}_{d.hour:02}-{d.minute:02}'.format(d=t0)

    file = f"{pth}{prefix}_{datestr}*.h5"

    file_list = glob.glob(file)[0]
    file_number = file_list.split('_')[-1]
    file_number = file_number.split('.')[0]
    file_number = int(file_number)
    return file_number

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
                      number_of_files=1):

    data = np.array([])
    time = np.array([])
    file_number = get_file_number(pth,file_base_name,t0)
    dt = datetime.timedelta(minutes=1) # Assume one minute file duration
    
    this_files_date = t0
    for i in range(number_of_files):

        this_file_number = file_number + i
#         date_str = this_files_date.strftime("%Y-%m-%d_%H-%M-%S")
        date_str = this_files_date.strftime("%Y-%m-%d_%H-%M") + "-00"
        filename = file_base_name + '_' + date_str + '_UTC_' + f'{this_file_number:06}' + '.h5'
        this_file = pth+filename
        try:
            f = h5py.File(this_file,'r')
            this_data = np.array(f['Acquisition/Raw[0]/RawData'][:,chan_min:chan_max])
            this_time = np.array(f['Acquisition/Raw[0]/RawDataTime'])
            
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
        
    return data, time, attrs