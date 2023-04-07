import numpy as np
import datetime
import h5py
import glob
from scipy.signal import detrend
from numpy.fft import fftshift, fft2, fftfreq
from datetime import datetime as DT
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import geopy.distance
from shapely.geometry import LineString
import pandas as pd
from scipy.signal import butter, filtfilt,detrend
from scipy.fft import fft, ifft
import scipy
import obspy
from obspy import UTCDateTime
from obspy.core import Trace
from obspy.clients.fdsn import Client
from libcomcat.search import search
from libcomcat.dataframes import get_summary_data_frame
from obspy.signal.trigger import classic_sta_lta,recursive_sta_lta
from obspy.signal.trigger import plot_trigger
from scipy.signal import butter, filtfilt,detrend
import numpy as np
from numpy.polynomial import legendre as leg
from numpy.linalg import lstsq
import datetime
from datetime import timedelta
import csv

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
    from obspy import UTCDateTime
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
                      verbose=False,
                      pad=False):

    data = np.array([])
    time = np.array([])

    
    dt = datetime.timedelta(minutes=1) # Assume one minute file duration
    this_files_date = t0
    
    for i in range(number_of_files):

        file_number = get_file_number(pth,file_base_name,this_files_date,verbose=verbose)
        if file_number == -1:
            raise ValueError('Failed to find file number.')
#             return [-1], [-1], [-1]
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
            #print('File problem with: %s'%this_file)
            #print(e)
            print('stuff')
            # There's probably a better way to handle this...
            #             return [-1], [-1], [-1]


        this_files_date = this_files_date + dt
    
    #if pad==True:
        # Add columns of zeros to give data matrix the correct dimensions
        
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
               extent=[x_lims[0],x_lims[-1],x_max,0])
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis_date()
    plt.grid()
    
    # Subplot: Single DAS Channel
    ax = plt.subplot(4,1,2)
    fig.patch.set_facecolor('w')
#     graph_spacing = -400
    graph_spacing = -20
    for jj in (41,400,800,1400):
        plt.plot(dates,datafilt[:,jj]-jj/graph_spacing,label=f'OD = {int(jj*dx)} m')
    plt.legend(loc='upper right')
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
    
    
def data_quicklook(     dates,datafilt,
                        x_max,stitle,filename=None,
                        das_vmax=0.1,
                        network_name='',
                        ylim=None):
    '''
    Make a nice plot of DAS data 
    '''
    dx = x_max / datafilt.shape[1]
    fig,ax=plt.subplots(figsize=(10,10))
    date_format = mdates.DateFormatter('%H:%M:%S')
    
    # Subplot: DAS Data

    ax.set_title(f'{network_name}')
    # plt.imshow(datafilt.T,vmin=-0.1,vmax=0.1,cmap='seismic',aspect='auto')
    x_lims = mdates.date2num(dates)
    plt.imshow(datafilt.T,vmin=-das_vmax,vmax=das_vmax,
               cmap='seismic',aspect='auto', 
               extent=[x_lims[0],x_lims[-1],x_max,0])
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis_date()
    plt.grid()
    if ylim is not None:
        ax.set_ylim(ylim)
    

    fig.suptitle(stitle,fontsize=20)
    plt.tight_layout()
    
    if filename==None:
        plt.show()
    else:
        plt.savefig(filename)
        plt.close()
        
def fk_analysis(t0, draw_figure = True,downsamplefactor=5,cable = 'whidbey', record_length = 1):
    
    prefix, network_name, datastore = data_wrangler(cable,record_length,t0)
    try:
        data,dates,attrs = open_sintela_file(prefix,
                                         t0,
                                         datastore,
                                         number_of_files=record_length,
                                         verbose=False)
    except:
        print("error'ed out")
        return [np.nan], [np.nan], [np.nan]
    
    x1 = 1225     # Entire subsea region
    x2 = 1600

    subsea_data = detrend(data[:,x1:x2])
    downsampled_subsea_data = subsea_data[::downsamplefactor,:]

    ft = fftshift(fft2(downsampled_subsea_data))
    f = fftshift(fftfreq(downsampled_subsea_data.shape[0], d=0.01 * downsamplefactor))
    k = fftshift(fftfreq(downsampled_subsea_data.shape[1], d=attrs['SpatialSamplingInterval']))
    
    return ft,f,k


#Function that returns the DAS data according to USGS earthquake events

def das_downloader(event_df, this_id, cab):
    this_event = event_df[event_df.id==this_id]
    t0 = this_event['time'].iloc[0]

    cable = cab
    record_length = 1 #minutes

    try:
        if cable == 'seadasn':
            prefix = 'seadasn'
            network_name = 'SeaDAS-N'
            if t0 < datetime.datetime(2022, 6, 20, 0, 0, 0):
                datastore='/data/data0/seadasn_2022-03-17_2022-06-20/'
            elif t0 > datetime.datetime(2022, 6, 20, 0, 0, 0) and t0 < datetime.datetime(2022, 10, 6, 0, 0, 0):
                datastore='/data/data7/seadasn_2022-06-21_2022-10-06/'
            elif t0 > datetime.datetime(2022, 11, 1, 0, 0, 0) and t0 < datetime.datetime(2022, 11, 5, 0, 0, 0):
                datastore='/data/data4/seadasn_2021-11-01_2021-11-05/'
            elif t0 > datetime.datetime(2022, 10, 13, 0, 0, 0) and t0 < datetime.datetime(2022, 11, 1, 0, 0, 0):
                datastore='/data/data1/seadasn_2021-10-13_2021-11-01/'

        elif cable == 'whidbey':
            prefix = 'whidbey'
            network_name='Whidbey-DAS'
            
            if t0 < datetime.datetime(2022, 10, 23, 4, 49, 0):
                datastore = '/data/data5/Converted/'
                
            elif t0 > datetime.datetime(2022, 10, 23, 4, 49, 0):
                datastore = '/data/data6/whidbey/'
                

        data,times,attrs = open_sintela_file(prefix,
                                             t0,
                                             datastore,
                                             number_of_files=record_length,
                                             verbose=True)
        x_max=data.shape[1] * attrs['SpatialSamplingInterval']
        #print('data:', data)
        low_cut = 3
        hi_cut = 8
        
        b,a = butter(2,[low_cut,hi_cut],'bp',fs=attrs['MaximumFrequency']*2,output='sos')
        data_filt = filtfilt(b,a,data,axis=0)

    except Exception as e:
        print(f'caught {type(e)}: e')
        data = None
        times = None
        dates = None
        attrs = None
        x_max = None
        data_filt = None
    return data, times, attrs, x_max, this_id, data_filt, t0


#For sintela fiberroute and calibration files
def fiber_channel_locator(das_data, attrs, fiber_file, cal_file, chan_spac = 6.3, num_chans = 1750, save_file = False):
    fiber_location = pd.read_csv(fiber_file, header=1)
    fiber_calibration = pd.read_csv(cal_file, header=1)
    fiber_distance = []
    opt_dis_merge = []
    chan_num = []
    

    
    for index_fib, row_fib in fiber_location.iterrows():
        if index_fib == 0 :
            fiber_distance.append(0)
        elif index_fib > 0:
            coords_1 = (row_fib['Lat'], row_fib['Long'])
            coords_2 = (fiber_location.iloc[index_fib-1]['Lat'], fiber_location.iloc[index_fib-1]['Long'])
            distance = geopy.distance.geodesic(coords_1, coords_2).m
            fiber_distance.append(distance   + fiber_distance[-1])

        l = []
        for index_cal, row_cal in fiber_calibration.iterrows():
            if row_fib['Lat'] == row_cal['Lat'] and row_fib['Long'] == row_cal['Long']:
                l.append(row_cal['Opt Dist'])
            else:
                pass
        if l:
            opt_dis_merge.append(l[0])
        else:
            opt_dis_merge.append(np.nan)

    fiber_location['Fiber Dist'] = fiber_distance
    fiber_location['Opt Dist'] = opt_dis_merge
    dis_interp = fiber_location['Opt Dist'].interpolate(method='linear', fill_value='extrapolate')
    fiber_location['Opt Dist Interp'] = dis_interp
    
    for index_merge, row_merge in fiber_location.iterrows():
        if row_merge['Opt Dist Interp'] != np.nan:
            
            chan_num.append(row_merge['Opt Dist Interp'] / chan_spac)
            
    fiber_location['Chan Num Interp'] = chan_num
    
    coords_of_chans_x = []
    coords_of_chans_y = []

    channel_number_counter = 0
    for index, values in fiber_location.iterrows():

        if index < len(fiber_location) - 1:
            
            xy_floats = [(fiber_location.iloc[index]['Long'], fiber_location.iloc[index]['Lat']),
                         (fiber_location.iloc[index+1]['Long'], fiber_location.iloc[index+1]['Lat'])]
            line = LineString(xy_floats)
            
            num_points = int(round((fiber_location.iloc[index+1]['Opt Dist Interp'] - values['Opt Dist Interp']) / attrs['SpatialSamplingInterval']))
            channel_number_counter += num_points
            #print(channel_number_counter)
            
            channel_difference = das_data.shape[1] - channel_number_counter
            
            if index == len(fiber_location) - 2:

                num_points = int(round((fiber_location.iloc[index+1]['Opt Dist Interp'] - values['Opt Dist Interp']) / attrs['SpatialSamplingInterval'])) + channel_difference
                
                new_points = [line.interpolate(i/float(num_points - 1), normalized=True) for i in range(num_points)]
                xs = [point.x for point in new_points]
                ys = [point.y for point in new_points]
                coords_of_chans_x.append(xs)
                coords_of_chans_y.append(ys)
            
            else:

                num_points = int(round((fiber_location.iloc[index+1]['Opt Dist Interp'] - values['Opt Dist Interp']) / attrs['SpatialSamplingInterval']))
                new_points = [line.interpolate(i/float(num_points - 1), normalized=True) for i in range(num_points)]
                xs = [point.x for point in new_points]
                ys = [point.y for point in new_points]
                coords_of_chans_x.append(xs)
                coords_of_chans_y.append(ys)
            


            
    

    flat_x =  [item for sublist in coords_of_chans_x for item in sublist] #Longitudes
    flat_y =  [item for sublist in coords_of_chans_y for item in sublist] #Latitudes
    
    if save_file == True:
        
        interp_channel_dict = {'Channel Number': np.arange(1, len(flat_x)+1, 1),
                              'Longitude': flat_x,
                              'Latitude': flat_y}
        
        chan_interp_frame = pd.DataFrame(data = interp_channel_dict)

        chan_interp_frame.to_csv('interpolated_channel_locations.csv', index=False)
        
    
    return fiber_location, fiber_calibration, flat_x, flat_y