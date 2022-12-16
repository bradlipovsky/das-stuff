import os
from tqdm import tqdm
from dasquakes import *
import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

# cable='whidbey'
cable = 'seadasn'
record_length = 10 
t0 = datetime.datetime(2022, 11, 11, 0, 16-1, 0)

# boioioing
lower_freq1 = .01
upper_freq1 = .1

lower_freq2 = .1
upper_freq2 = 1

# boioioing channels
lower_chan = 65
upper_chan = 80

vm = 50

counter = 0
for day in range(11, 15, 1):
    for hour in tqdm(range(0, 24, 1)):
    #     if counter == 1:
    #         break
        for minute in range(0, 60, record_length):
            t0 = datetime.datetime(2022, 11, day, 0 + hour, minute, 0)
    #         print(t0)


            prefix, network_name, datastore = data_wrangler(cable, record_length, t0)

            data, dates, attrs = open_sintela_file(
                 prefix,
                 t0,
                 datastore,
                 number_of_files=record_length,
                 verbose=False
            )

            nt = data.shape[0]
            nx = data.shape[1]
            data_subset = data[:, lower_chan:upper_chan]
            fs = 2 * attrs['MaximumFrequency']


            b,a = butter(2,[lower_freq1, upper_freq1],'bandpass',fs = fs)
            d,c = butter(2,[lower_freq2, upper_freq2],'bandpass',fs = fs)

            data_subset_filtered1 = filtfilt(b,a,data_subset,axis=0)
            data_subset_filtered2 = filtfilt(d,c,data_subset,axis=0)
            y_max = len(data_subset) / 100



            time_change = datetime.timedelta(hours = -8)
            minute_change = datetime.timedelta(minutes = record_length)
            new_time = t0 + time_change
            time_end = t0 + time_change + minute_change



            if data_subset.max() >= 30:
                fig, axs = plt.subplots(
                    2, 3, figsize = (18,18), gridspec_kw={'height_ratios': [3, 1]}
                )
                # unfiltered data
                axs[0,0].imshow(
                    data_subset,
                    aspect = 'auto',
                    vmin = -vm,
                    vmax = vm,
                    extent = (lower_chan, upper_chan, y_max, 0)
                )
                axs[0, 0].set_title(
                        'Raw Data'  '\n' +
                        str(new_time) +
                        ' - ' + 
                        str(time_end.time()),
                        fontsize = '15'
                    )
                axs[0,0].tick_params(axis='both', labelsize=15)
                axs[0,0].set_ylabel('Timestep (1/100 seconds)', fontsize = '15')
                axs[0,0].set_xlabel('Channel', fontsize = '15')


                # filtered data
                axs[0,1].imshow(
                    data_subset_filtered1,
                    aspect = 'auto',
                    vmin = -vm,
                    vmax = vm,
                    extent = (lower_chan, upper_chan, y_max, 0)
                )
                axs[0, 1].set_title(
                        'Filtered Data ' + str(lower_freq1) + ' - ' + 
                        str(upper_freq1) + ' Hz\n' +
                        str(new_time) +
                        ' - ' + 
                        str(time_end.time()),
                        fontsize = '15'
                    )
                axs[0, 1].set_ylabel('Timestep (1/100 seconds)', fontsize = '15')
                axs[0, 1].set_xlabel('Channel', fontsize = '15')
                axs[0,1].tick_params(axis='both', labelsize=15)


                # filtered data
                axs[0,2].imshow(
                    data_subset_filtered2,
                    aspect = 'auto',
                    vmin = -vm,
                    vmax = vm,
                    extent = (lower_chan, upper_chan, y_max, 0)
                )
                axs[0, 2].set_title(
                        'Filtered Data ' + str(lower_freq2) + ' - ' + 
                        str(upper_freq2) + ' Hz\n' +
                        str(new_time) +
                        ' - ' + 
                        str(time_end.time()),
                        fontsize = '15'
                    )
                axs[0, 2].set_ylabel('Timestep (1/100 seconds)', fontsize = '15')
                axs[0, 2].set_xlabel('Channel', fontsize = '15')
                axs[0,2].tick_params(axis='both', labelsize=15)


                timestamp_range = (pd.date_range(start=new_time,
                                                end = time_end,
                                               freq = 'min'))
                time_labels = []

                for i in timestamp_range:
                    time_labels.append(str(i)[11:-3])

                N = 11
                ind = np.arange(N)      


                # unfiltered time series
                axs[1,0].plot(dates ,data_subset[:,9])

                axs[1,0].set_title(
                        'Channel 74 Time Series',
                        fontsize = '15'
                    )
                axs[1,0].tick_params(axis='both', labelsize=15)
                axs[1,0].set_xticklabels(time_labels)
                axs[1,0].tick_params(axis='x', labelrotation = 45)


                # filtered time series
                axs[1,1].plot(dates, data_subset_filtered1[:,9])
                axs[1,1].set_title('Channel 74 Time Series', fontsize = '15')
    #             axs[1,1].tick_params(axis='both', labelsize=15)
                axs[1,1].tick_params(axis = 'both', labelsize = 15)
                axs[1,1].set_xticklabels(time_labels)
                axs[1,1].tick_params(axis='x', labelrotation = 45)


                            # filtered time series
                axs[1,2].plot(dates, data_subset_filtered2[:,9])
                axs[1,2].set_title('Channel 74 Time Series', fontsize = '15')
    #             axs[1,1].tick_params(axis='both', labelsize=15)
                axs[1,2].tick_params(axis = 'both', labelsize = 15)
                axs[1,2].set_xticklabels(time_labels)
                axs[1,2].tick_params(axis='x', labelrotation = 45)

                fig.patch.set_facecolor('w')

                plt.locator_params(axis='y', nbins=15)
                axs[0,0].grid(color = 'red', which = 'both')
                axs[0,1].grid(color = 'red', which = 'both')
                axs[0,2].grid(color = 'red', which = 'both')

                plt.tight_layout()
                plt.savefig('figs/' + str(new_time) + ' -> ' + str(time_end) + '.png')
#                 plt.savefig('figs/' + str(new_time) + ' -> ' + str(time_end) + '.eps')