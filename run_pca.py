'''
SVD analysis of Distributed Acoustic Sensing (DAS) frequency-wavenumber (FK) plots.

TODO: Load acquisition parameters from each file rather than assuming they are fixed for the entire duration.

TODO: Add wind speed data

TODO: Convert frequency-wavenumber to frequency-wavespeed

TODO: Experiment with stacking parameters.  Note that SGWs appear with longer stacking windows.

TODO: Speed up the calculations using multiprocessing.
'''
import h5py
import matplotlib.pyplot as plt
import numpy as np
import datetime
from time import perf_counter
from scipy.sparse.linalg import svds
from tqdm import tqdm
from dasquakes import *
import pickle

def main():
    t0 = perf_counter()
    '''
    Parameters for the analysis
    '''

    q = 10           # decimation factor
    N = 24*100       # number of steps
    dt = 60          # minutes between steps
    record_length=15 # minutes per step
    nt = int(record_length*6000/q) # Number of time steps in each sample
    nx = 375         # Number of subsea channels at Whidbey
    filename = 'svd.pickle'
    
    '''
    Begin the workflow
    '''
    svd_analysis(N=N,dt=dt,q=q,outputfile=filename,record_length=record_length)
    
    file = open(filename, 'rb')
    U,S,V,t,f,k = pickle.load(file)
    file.close()
    
    f = fftshift(fftfreq(nt, d=0.01 * q))
    k = fftshift(fftfreq(nx, d=6.38))
    
    first_mode = U[:,5]
    first_mode = np.abs(U[:,5].reshape((nt,nx))/np.max(np.abs(first_mode)))
    first_time_series = V[5,:]
    
    second_mode = U[:,4]
    second_mode = np.abs(U[:,4].reshape((nt,nx))/np.max(np.abs(second_mode)))
    second_time_series = V[4,:]
    
    plot_svd(S,f,k,t,first_mode,first_time_series,'svd_plot_mode-1.pdf')
    plot_svd(S,f,k,t,second_mode,second_time_series,'svd_plot_mode-1.pdf')
    print(f'Total runtime: {perf_counter()-t0} s')
    
    
def svd_analysis(q=10,N=24,dt=60,record_length=2,
                 start_time = datetime.datetime(2022, 5, 8, 0, 0, 0), 
                 outputfile='svd.pickle'):
    
    
    '''
    Build the data matrix
    '''
    fs = 100
    file_duration = 60
    nt = int(record_length*file_duration*fs/q) # Number of time steps in each sample
    nx = 375 # Number of subsea channels at Whidbey
    D = np.zeros((nx*nt,N))
    t = []
    
    for i in tqdm(range(N)):
        this_time = start_time + i*datetime.timedelta(minutes=dt)
        t.append(this_time)
        ft,f,k = fk_analysis(this_time,draw_figure=False,downsamplefactor=q,
                            record_length = record_length)
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
    file = open(outputfile, 'wb')
    pickle.dump((U,S,V,t,f,k), file)
    file.close()
    
    

def plot_svd(S,f,k,t,mode,time_series,filename='svd_plot.pdf'):

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
    
    plt.savefig(filename)


if __name__ == "__main__":
    main()