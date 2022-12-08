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


def dt_to_utc_format(t):
    from obspy import UTCDateTime
    return UTCDateTime(t.strftime('%Y-%m-%dT%H:%M:%S'))

def utc_to_dt_format(t):
    dt_str = t.strftime('%Y/%m/%d %H:%M:%S')
    format1  = "%Y/%m/%d %H:%M:%S"
    dt_utc = DT.strptime(dt_str, format1)
    return dt_utc
    
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
        
def fk_analysis(t0, 
                draw_figure = True,
                fs=10, dx=6.38,
                cable = 'whidbey', 
                record_length = 1,
                distance_range=[7816,10208],
                verbose = True):
    '''
    This function takes inputs that describe a subset of a DAS deployment and returns FK data.

    The default distance range represents the subsea part of the whidbey cable
    '''
    
    prefix, network_name, datastore = data_wrangler(cable,record_length,t0)
    try:
        data,dates,attrs = open_sintela_file(prefix,
                                         t0,
                                         datastore,
                                         number_of_files=record_length,
                                         verbose=True)
    except:
        print("error'ed out")
        return [np.nan], [np.nan], [np.nan]
    
    ''' 
    Downsampling
    '''
    fs_input = 2*attrs['MaximumFrequency']
    t_downsample_factor = int(fs_input/fs)
    if verbose: print(f"Temporal downsampling factor = {t_downsample_factor}")
    
    dx_input = attrs['SpatialSamplingInterval']
    x_downsample_factor = int(dx_input/dx)
    if verbose: print(f"Spatial downsampling factor = {x_downsample_factor}")
    
    x1 = int(distance_range[0]/dx_input)
    x2 = int(distance_range[1]/dx_input)

    subsea_data = detrend(data[:,x1:x2:x_downsample_factor])
    downsampled_subsea_data = subsea_data[::t_downsample_factor,:]
    if verbose: print(f"Data dimension: {downsampled_subsea_data.shape}")
    
    '''
    Calculate FFT
    '''

    ft = fftshift(fft2(downsampled_subsea_data))
    f = fftshift(fftfreq(downsampled_subsea_data.shape[0], d=1/fs))
    k = fftshift(fftfreq(downsampled_subsea_data.shape[1], d=dx))
    
    return ft,f,k

def plot_svd(S,f,k,t,mode,time_series,outputfile='svd'):

    '''
    Plot the results
    '''
    vm = 0.1
    
    plt.subplots(2,1,figsize=(10,10))

    ax1=plt.subplot(2,1,1)
    plt.title(f'Fraction of variance in 1st mode: {100*max(S)/sum(S)}%')
    c=plt.imshow(mode,aspect='auto',vmin=0,vmax=vm,extent=[k[0],k[-1],f[0],f[-1]],cmap='gray_r')

    ax1.set_ylim([-2.5,2.5])
    ax1.set_xlim([-0.04,0.04])
    ax1.set_xlabel('Wavenumber (1/m)')
    ax1.set_ylabel('Frequency (Hz)')
#     plt.colorbar()

    ax2=plt.subplot(2,1,2)
    ind = np.where(np.abs(time_series)>1e-10)
    sign_change = np.sign(np.mean(time_series))
    ax2.plot(t[ind],time_series[ind]*sign_change,'o')
    plt.xticks(rotation = 25)
    ax2.grid()
    
    plt.savefig(f'{outputfile}.pdf')


def svd_analysis(q=10,N=24,dt=60,
                 start_time = datetime.datetime(2022, 5, 8, 0, 0, 0), 
                 outputfile='svd_10min_window_length',
                 window_length=1,
                 verbose=False):
    '''
    Carries out an SVD analysis of a data matrix, D. Each column of D contains a flattened fk-plot.  There are N
    samples contained in D.  The samples are spaced apart in time by dt minutes.  For each sample, we evaluate chunks of
    window_length minutes of data. q is the decimation factor.
    '''
    
    '''
    Build the data matrix
    '''
    sample_rate = 100 # Hz
    file_duration = 60 # s
    samples_per_minute = sample_rate * file_duration
    nt = int(window_length * samples_per_minute/q) # Number of time steps in each sample
    nx = 375 # Number of subsea channels at Whidbey
    D = np.zeros((nx*nt,N))
    t = []
    
    if verbose: print(f'Building data matrix with nt={nt}, nx={nx}') 
    for i in tqdm(range(N)):
        this_time = start_time + i*datetime.timedelta(minutes=dt)
        t.append(this_time)
        ft,f,k = fk_analysis(this_time,draw_figure=False,downsamplefactor=q,
                            record_length = window_length)
        if len(ft) == 1:
            continue

        shape = ft.shape
        this_nt = shape[0]
        this_nx = shape[1]

        if this_nt < nt:
            ft_new = np.zeros((nt,nx))
            ft_new[0:this_nt,0:nx] = np.abs(ft)
            this_column =  ft_new.flatten()
        elif this_nt > nt:
            ft_new = np.zeros((nt,nx))
            ft_new[0:nt,0:nx] = np.abs(ft[0:nt,0:nx])
            this_column =  ft_new.flatten()
        else:
            this_column = np.abs( ft.flatten() )

        D[:,i] = this_column
    t=np.array(t)

    '''
    Calculate the SVD
    '''
    ns = N
    t1 = perf_counter()
    U,S,V = svds( D[:,0:ns] )
    t = t[0:ns]
    print(f'SVD runtime:   {perf_counter()-t1} s')
    
    # open a file, where you ant to store the data
    file = open(f'{outputfile}.pickle', 'wb')
    pickle.dump((U,S,V,t,f,k), file)
    file.close()

    
    
def svd_analysis(N=24,dt=60,dx=6.38,fs=10,
                 distance_range=[7816,10208], 
                 record_length=2,
                 start_time = datetime.datetime(2022, 5, 8, 0, 0, 0), 
                 outputfile='svd.pickle',
                 verbose=False):
    '''
    Build the data matrix
    '''
    file_duration = 60 #seconds.  this shouldn't change through the deployment.

    # Number of time steps in each sample
    nt = int(record_length*file_duration*fs) 
    
    # Number of subsea channels at Whidbey. This also changes.
    nx = int((distance_range[1]-distance_range[0]+1)/dx) 
    
    if verbose: print(f"nx={nx}, nt={nt}")
    
    D = np.zeros((nx*nt,N))
    t = []

    for i in tqdm(range(N)):

        this_time = start_time + i*datetime.timedelta(minutes=dt)
        t.append(this_time)
        ft,f,k = fk_analysis(this_time,draw_figure=False,fs=fs, 
                            record_length = record_length,
                            distance_range = distance_range,
                            verbose=verbose)
        if len(ft) == 1:
            continue

        this_nt = ft.shape[0]
        this_nx = ft.shape[1]

        if this_nt < nt:
            ft_new = np.zeros((nt,nx))
            ft_new[0:this_nt,0:nx] = np.abs(ft)
            this_column =  ft_new.flatten()
        elif this_nt > nt:
            ft_new = np.zeros((nt,nx))
            ft_new[0:nt,0:nx] = np.abs(ft[0:nt,0:nx])
            this_column =  ft_new.flatten()
        else:
            # This should be the typically situation: there are exactly as many time steps as we anticipated.
            this_column = np.abs( ft.flatten() )

        D[:,i] = this_column


    t=np.array(t)
    
    '''
    Calculate the SVD
    '''
    ns = N
    t1 = perf_counter()
    U,S,V = svds( D[:,0:ns] )
    t = t[0:ns]
    print(f'SVD runtime:   {perf_counter()-t1} s')
    
    # open a file, where you ant to store the data
    file = open(outputfile, 'wb')
    pickle.dump((U,S,V,t,f,k,nt,nx), file)
    file.close()
    
    

def plot_svd(f,k,t,mode,time_series,var,filename='svd_plot.pdf'):

    '''
    Plot the results
    '''
    vm = 0.1
    
    plt.subplots(2,1,figsize=(10,10))

    ax1=plt.subplot(2,1,1)
    plt.title(f'Fraction of variance in this mode: {var}%')
    c=plt.imshow(mode,aspect='auto',vmin=0,vmax=vm,extent=[k[0],k[-1],f[0],f[-1]],cmap='gray_r')

    ax1.set_ylim([-2.5,2.5])
    ax1.set_xlim([-0.04,0.04])
    ax1.set_xlabel('Wavenumber (1/m)')
    ax1.set_ylabel('Frequency (Hz)')
#     plt.colorbar()

    ax2=plt.subplot(2,1,2)
    ind = np.where(np.abs(time_series)>1e-10)
    sign_change = np.sign(np.mean(time_series))
    ax2.plot(t[ind],time_series[ind]*sign_change,'o')
    plt.xticks(rotation = 25)
    ax2.grid()
    
    plt.savefig(filename)