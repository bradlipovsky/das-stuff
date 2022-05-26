import numpy as np
import datetime
import h5py

import glob
def get_file_number(pth,t0):
    
    datestr = '{d.year}-{d.month:02}-{d.day:02}_{d.hour:02}-{d.minute:02}'.format(d=t0)

    file = f"{pth}seadasn_{datestr}*.h5"

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

def open_sintela_file(file_base_name,number_of_files,sample_rate,
                     file_duration,t0,file_number,pth,dt,
                             chmin,chanmax):

    data = np.array([])
    time = np.array([])
    samples_per_file = round(sample_rate * file_duration * 60)

    this_files_date = t0
    for i in range(number_of_files):

        this_file_number = file_number + i
        date_str = this_files_date.strftime("%Y-%m-%d_%H-%M-%S")
        filename = file_base_name + '_' + date_str + '_UTC_' + f'{this_file_number:06}' + '.h5'
        this_file = pth+filename
        try:
            f = h5py.File(this_file,'r')
            if i == 0:
                data = f['Acquisition/Raw[0]/RawData'][:,chmin:chanmax] 
                time = f['Acquisition/Raw[0]/RawDataTime']
            else:
                data = np.concatenate((data, f['Acquisition/Raw[0]/RawData'][:,chmin:chanmax] ))
                time = np.concatenate((time,f['Acquisition/Raw[0]/RawDataTime']))
        except Exception as e: 
            print('File problem with: %s'%this_file)
            print(e)


        this_files_date = this_files_date + dt
        
    return data, time