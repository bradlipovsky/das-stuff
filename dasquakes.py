import numpy as np
import datetime
import h5py

import glob

import obspy
from obspy import UTCDateTime
from datetime import datetime as DT

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def data_wrangler(cable,record_length,t0):
    if cable == 'seadasn':
        prefix = 'seadasn'
        network_name = 'SeaDAS-N'
        if t0 < datetime.datetime(2022, 6, 20, 0, 0, 0):
            datastore='/data/data0/seadasn_2022-03-17_2022-06-20/'
        else:
            datastore='/data/data7/seadasn/'

    elif cable == 'whidbey':
        prefix = 'whidbey'
        network_name='Whidbey-DAS'
        datastore = '/data/data5/Converted/'
        
    return prefix, network_name, datastore


def dt_to_utc_format(t):
    return UTCDateTime(t.strftime('%Y-%m-%dT%H:%M:%S'))

def utc_to_dt_format(t):
    dt_str = t.strftime('%Y/%m/%d %H:%M:%S')
    format1  = "%Y/%m/%d %H:%M:%S"
    dt_utc = DT.strptime(dt_str, format1)
    return dt_utc


def get_file_number(pth,prefix,t0,verbose=False):
    
    datestr = '{d.year}-{d.month:02}-{d.day:02}_{d.hour:02}-{d.minute:02}'.format(d=t0)

    file = f"{pth}{prefix}_{datestr}*.h5"
    if verbose:
        print(file)

    if len(glob.glob(file)) > 0:
        file_list = glob.glob(file)[0]
#         print(glob.glob(file))
        file_number = file_list.split('_')[-1]
        file_number = file_number.split('.')[0]
        file_number = int(file_number)
        return file_number
    else:
        return -1
    
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
                      verbose=False):

    data = np.array([])
    time = np.array([])

    
    dt = datetime.timedelta(minutes=1) # Assume one minute file duration
    this_files_date = t0
    
    for i in range(number_of_files):

        file_number = get_file_number(pth,file_base_name,this_files_date,verbose=verbose)
        if file_number == -1:
            print('Failed to find file number.')
            return [-1], [-1], [-1]
        date_str = this_files_date.strftime("%Y-%m-%d_%H-%M") + "-00"
        this_file = f'{pth}{file_base_name}_{date_str}_UTC_{file_number:06}.h5'
        
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
            
            # There's probably a better way to handle this...
#             return [-1], [-1], [-1]


        this_files_date = this_files_date + dt
        
    return data, time, attrs

def local_earthquake_quicklook(dates,datafilt,st,st2,
                        x_max,stitle,filename=None,
                        skip_seismograms=False,
                        das_vmax=0.1,
                        network_name=''):
    '''
    Make a nice plot of the DAS data and some local seismic stations
    '''
    dx = x_max / datafilt.shape[1]
    fig,ax=plt.subplots(figsize=(8,12))
    date_format = mdates.DateFormatter('%H:%M:%S')
    
    # Subplot: DAS Data
    ax=plt.subplot(4,1,1)
    ax.set_title(f'{network_name}')
    # plt.imshow(datafilt.T,vmin=-0.1,vmax=0.1,cmap='seismic',aspect='auto')
    x_lims = mdates.date2num(dates)
    plt.imshow(datafilt.T,vmin=-das_vmax,vmax=das_vmax,
               cmap='seismic',aspect='auto', 
               extent=[x_lims[0],x_lims[-1],0,x_max])
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis_date()
    plt.grid()
    
    # Subplot: Single DAS Channel
    ax = plt.subplot(4,1,2)
    fig.patch.set_facecolor('w')
    graph_spacing = -400
    for jj in (41,400,800,1400):
        plt.plot(dates,datafilt[:,jj]-jj/graph_spacing,label=f'OD = {int(jj*dx)} m')
    plt.legend(loc='upper left')
    ax.set_title(f'{network_name} Individual Channels')
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis_date()
    ax.autoscale(enable=True, axis='x', tight=True)
    plt.grid()


    if skip_seismograms==False:
        
        # Subplot:  station 1
        ax = plt.subplot(4,1,3)
        for tr in st:
            times_from_das = np.linspace(x_lims[0],x_lims[-1],len(tr.data))
            plt.plot(times_from_das,tr.data)
        fig.patch.set_facecolor('w')
        ax.set_title('UW NOWS HNN')
        ax.xaxis.set_major_formatter(date_format)
        ax.xaxis_date()
        ax.set_xlim((min(times_from_das),max(times_from_das)))
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
        ax.set_xlim((min(times_from_das),max(times_from_das)))
        plt.grid()
    
    

    fig.suptitle(stitle,fontsize=20)
    plt.tight_layout()
    
    if filename==None:
        plt.show()
    else:
        plt.savefig(filename)
        plt.close()
    
    
