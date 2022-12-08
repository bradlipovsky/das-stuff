'''
SVD analysis of Distributed Acoustic Sensing (DAS) frequency-wavenumber (FK) plots.

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
    Carry out the SVD analysis
    
    default parameters:
    svd_analysis(N=24,dt=60,dx=6.38,fs=10,
                     distance_range=[7816,10208], 
                     record_length=2,
                     start_time = datetime.datetime(2022, 5, 8, 0, 0, 0), 
                     outputfile='svd.pickle',
                     verbose=False):
    '''
    filename = 'svd.pickle'
    svd_analysis(record_length=5,N=24*262,
                start_time = datetime.datetime(2022, 3, 22, 0, 0, 0))
    
    '''
    Load the result and make plots
    '''
    file = open(filename, 'rb')
    U,S,V,t,f,k,nt,nx = pickle.load(file)
    file.close()
    
    for i in range(2):
        normalization = np.max(np.abs(U[:,5-i]))
        mode = np.abs(U[:,5-i].reshape((nt,nx))) / normalization
        time_series = V[5-i,:]
        variance = 100*S[5-i]/sum(S)
        plot_svd(f,k,t,mode,time_series,variance,f'svd_plot_mode-{i+1}.pdf')
    
    print(f'Total runtime: {perf_counter()-t0} s')

if __name__ == "__main__":
    main()