import numpy as np
import datetime
import h5py

import glob

import obspy
from obspy import UTCDateTime
from datetime import datetime as DT

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def dt_to_utc_format(t):
    return UTCDateTime(t.strftime('%Y-%m-%dT%H:%M:%S'))

def utc_to_dt_format(t):
    dt_str = t.strftime('%Y/%m/%d %H:%M:%S')
    format1  = "%Y/%m/%d %H:%M:%S"
    dt_utc = DT.strptime(dt_str, format1)
    return dt_utc


def get_file_number(pth,prefix,t0):
    
    datestr = '{d.year}-{d.month:02}-{d.day:02}_{d.hour:02}-{d.minute:02}'.format(d=t0)

    file = f"{pth}{prefix}_{datestr}*.h5"

    print(file)
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

def local_earthquake_quicklook(dates,datafilt,st,st2,x_max,event_df,catalog_index):

    
    
    fig,ax=plt.subplots(figsize=(8,12))
    date_format = mdates.DateFormatter('%H:%M:%S')
    
    # Subplot: DAS Data
    ax=plt.subplot(4,1,1)
    ax.set_title('SeaDAS-N')
    # plt.imshow(datafilt.T,vmin=-0.1,vmax=0.1,cmap='seismic',aspect='auto')
    x_lims = mdates.date2num(dates)
    plt.imshow(datafilt.T,vmin=-.1,vmax=.1,cmap='seismic',aspect='auto', extent=[x_lims[0],x_lims[-1],0,x_max])
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis_date()
    plt.grid()
    
    # Subplot: Single DAS Channel
    ax = plt.subplot(4,1,2)
    fig.patch.set_facecolor('w')
    plt.plot(dates,datafilt[:,800])
    ax.set_title('SeaDAS-N One Channel')
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis_date()
    plt.grid()

    
    
    # Subplot:  station 1
    ax = plt.subplot(4,1,3)
    for tr in st:
        times_from_das = np.linspace(x_lims[0],x_lims[-1],len(tr.data))
        plt.plot(times_from_das,tr.data)
    fig.patch.set_facecolor('w')
    ax.set_title('UW NOWS HNN')
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis_date()
    plt.grid()
    

    # Subplot:  station 2
    ax = plt.subplot(4,1,4)
    for tr in st2:
        times_from_das = np.linspace(x_lims[0],x_lims[-1],len(tr.data))
        plt.plot(times_from_das,tr.data)
    fig.patch.set_facecolor('w')
    ax.set_title('IU COR BH1')
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis_date()
    plt.grid()
    
    stitle=f"M{event_df.iloc[catalog_index]['Magnitude']}, {event_df.iloc[catalog_index]['Time UTC']} UTC"

    fig.suptitle(stitle,fontsize=20)
    plt.tight_layout()
    plt.show()