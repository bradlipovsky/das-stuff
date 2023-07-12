import h5py
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
import numpy as np
from numpy.fft import fftshift, fft2, fftfreq
import scipy


def correlate_das(s1,s2,mode="channel"):

    # throw an error of input sizes are inconsistent
    if s1.shape != s2.shape:
        raise ValueError("s1 and s2 must have the same size!")
    
    # get fft size
    sz = s1.shape[0]
    n_bits = 1+int(np.log2(2*sz-1))
    fft_sz = 2**n_bits
    
    # take FFT along time axis for both
    fft_s1 = np.fft.fft(s1, fft_sz, axis=0)
    fft_s2 = np.fft.fft(s2, fft_sz, axis=0)
        
    # take complex conjugate of second signal 
    fft_s2_conj = np.conj(fft_s2)
      
    # multiply to get correlation function
    corr_fft = fft_s1*fft_s2_conj
    
    # take inverse fourier transform
    corr = np.fft.ifft(corr_fft, axis=0)
    
    # normalize (this is taken from tslearn, and it is confusing). "Overall" option gives a single correlation 
    # function that represents all channels; "channel" option returns a separate correlation function 
    # for each channel of input data
    if mode == "overall":
        norm1 = np.linalg.norm(s1)
        norm2 = np.linalg.norm(s2)
        norm_factor = norm1*norm2
        corr = np.vstack((corr[-(sz-1) :], corr[:sz]))
        norm_corr = np.real(corr).sum(axis=-1) / norm_factor
        
    if mode == "channel":
        norm1 = np.linalg.norm(s1,axis=0)
        norm2 = np.linalg.norm(s2,axis=0)
        norm_factor = norm1*norm2
        corr = np.vstack((corr[-(sz-1) :], corr[:sz]))
        norm_corr = np.real(corr) / norm_factor

    # throw an error if an invalid mode is requested
    if mode != "overall" and mode != "channel":
        raise ValueError("Valid modes are channel and overall.")

    return norm_corr


def template_match_das(template,data,window_length): #xcorr_thresh):

    # define container
    all_corr = []

    num_windows = int(data.shape[0]/window_length)

    # iterate through time windows
    for i in range(num_windows):

        # pull out a time window of data
        # Define the start and end indices of the time window
        start_index = i*window_length
        end_index = start_index + window_length
        timewindow = data[start_index:end_index,:]
        
        # call cross correlation function
        corr = correlate_das(template,timewindow,mode="overall")
        #corr = scipy.signal.correlate(template,timewindow)

        # get maximum value
        max_corr = np.max(np.abs(corr))
        
        
        # save value
        all_corr.append(max_corr)
        
    return all_corr
    