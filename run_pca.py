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

def main():
    t0 = perf_counter()

    '''
    Parameters for the analysis
    '''
    start_time = datetime.datetime(2022, 5, 8, 0, 0, 0)
    q = 10  # decimation factor
    N = 24*100 # number of samples to analyze
    dt = 60 # number of minutes between samples
    
    '''
    Build the data matrix
    '''
    nt = int(6000/q) # Number of time steps in each sample
    nx = 375 # Number of subsea channels at Whidbey
    D = np.zeros((nx*nt,N))
    t = []
    
    for i in tqdm(range(N)):
        this_time = start_time + i*datetime.timedelta(minutes=dt)
        t.append(this_time)
        ft,f,k = fk_analysis(this_time,draw_figure=False,downsamplefactor=q)
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

    vm = 0.1
    first_mode = U[:,5]
    first_mode = np.abs(U[:,5].reshape((nt,nx))/np.max(np.abs(first_mode)))
    first_time_series = V[5,:]
    
    '''
    Plot the results
    '''
    plt.subplots(2,1,figsize=(10,10))

    ax1=plt.subplot(2,1,1)
    plt.title(f'Fraction of variance in 1st mode: {100*max(S)/sum(S):.1f}%')
    c=plt.imshow(first_mode,aspect='auto',vmin=0,vmax=vm,extent=[k[0],k[-1],f[0],f[-1]],cmap='gray_r')

    ax1.set_ylim([-2.5,2.5])
    ax1.set_xlim([-0.04,0.04])
    ax1.set_xlabel('Wavenumber (1/m)')
    ax1.set_ylabel('Frequency (Hz)')
#     plt.colorbar()

    ax2=plt.subplot(2,1,2)
    print(first_time_series)
    ind = np.where(np.abs(first_time_series)>1e-10)
    sign_change = np.sign(np.mean(first_time_series))
    ax2.plot(t[ind],first_time_series[ind]*sign_change,'o')
    plt.xticks(rotation = 25)
    ax2.grid()
    
    plt.savefig('pca.pdf')


    print(f'Total runtime: {perf_counter()-t0} s')


if __name__ == "__main__":
    main()