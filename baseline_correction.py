#### Script made for baseline correction in Seismosignal ######

import numpy as np
import scipy as sp
from scipy.signal import butter, lfilter
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt 
import matplotlib
from scipy import integrate
from scipy.fftpack import fft


def butter_bandpass(lowcut, highcut, fs, order=4):  ##  4th order Butterworth filter bandpass
    
    nyq = 0.5 * fs
    
    low = lowcut / nyq
    
    high = highcut / nyq
    
    b, a = butter(order, [low, high], btype='band')
    
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=4):  

    b, a = butter_bandpass(lowcut, highcut, fs, order=order)

    y = lfilter(b, a, data)

    return y


def time_series_to_FFT(time_series_var, time):

	No_time_step = len(time); 

	dt = time[1] - time[0];  

	if No_time_step % 2 == 0:
		
		print "Number of time steps is even number"; 

		N = int(No_time_step/2); 

		f= sp.arange(0, N+1, 1)*1.0/N*(0.5/dt);

	else: 

		print "Number of time steps is odd number"; 

		N = int((No_time_step-1)/2);

		f = sp.arange(0, N+1, 1)*1.0/(2*N+1)*(1/dt);  

	fft_var_full = fft(time_series_var); 

	fft_var = fft_var_full[0:(N+1)];  ### only half is the effective info 

	fft_var_mag = 2*abs(fft_var)/No_time_step; 

	return f, fft_var, fft_var_mag; 

def baseline_correction_filtering(data_x, data_y, lowcut, highcut, fs, order, Is_detrend, Is_bandpass_filter):   ## perform linear detrend and bandpass filtering for acc data in order to baseline correct

	if Is_detrend: 

		coef = np.polyfit(data_x, data_y, 1);  ### linear fit 

		poly1d_fn = np.poly1d(coef) ### linear fit function 

		data_y_detrend = data_y - poly1d_fn(data_x); 

	else: 
	
		data_y_detrend = data_y; 

	if Is_bandpass_filter: 

		data_y_filtered_detrend = butter_bandpass_filter(data_y_detrend, lowcut, highcut, fs, order);

	else: 

		data_y_filtered_detrend = data_y_detrend; 

	return data_y_filtered_detrend;  



def acc_to_dis(acc_x, acc_y):   ### double integration from acc to vel and dis, acc_x is time vector, acc_y is acceleration series  

	vel = integrate.cumtrapz(acc_y, acc_x, initial=0);  

	dis = integrate.cumtrapz(vel, acc_x, initial=0); 

	return vel, dis; 



if __name__ == '__main__':
	
	#### simple test case #####

	input_acc_file = 'wavefield_1D_x_from_depth_0.0_at_depth_12_acc.txt';  ## 'wavefield_1D_x_from_depth_0.0_at_depth_12_acc.txt'; 

	acc = np.loadtxt(input_acc_file);   ### two column, first column is time, second column is acceleration  

	dt = 0.02; 

	fs = 1/dt; 

	f_lowcut = 0.1; 

	f_highcut = 24.999;   ### Important Note:: To make filter take effect, f2 must be less than half of the samping frequency. i.e., Nyquist frequency 

	order = 4;

	Is_detrend = True; 

	Is_bandpass_filter = True; 

	#######################################################################################
	#######################################################################################
	################################ old script for backup purpose ########################
	#######################################################################################
	#######################################################################################

	# coef = np.polyfit(acc[:, 0], acc[:, 1], 1); 

	# # print coef; 

	# # coef = [-1e-5, 1.4e-4]; 

	# poly1d_fn = np.poly1d(coef) ### linear fit function 

	# acc_detrend = acc[:, 1] - poly1d_fn(acc[:, 0]); 

	# acc_filtered_detrend = butter_bandpass_filter(acc_detrend, f_lowcut, f_highcut, fs, order=4); 

	# # acc_filtered_detrend = butter_bandpass_filter(acc[:, 1], f_lowcut, f_highcut, fs, order=4); 

	# # Vel = integrate.cumtrapz(acc[:, 1], acc[:, 0], initial=0);  

	# # Vel = integrate.cumtrapz(acc_detrend, acc[:, 0], initial=0);  

	# Vel = integrate.cumtrapz(acc_filtered_detrend, acc[:, 0], initial=0);  

	# Dis = integrate.cumtrapz(Vel, acc[:, 0], initial=0); 

	# f_original, fft_var_original, fft_var_mag_original = time_series_to_FFT(acc[:, 1], acc[:, 0]);

	# f_filtered, fft_var_filtered, fft_var_mag_filtered = time_series_to_FFT(acc_filtered_detrend, acc[:, 0]);

	#######################################################################################
	#######################################################################################
	#######################################################################################
	#######################################################################################


	acc_filtered_detrend = baseline_correction_filtering(acc[:, 0], acc[:, 1], f_lowcut, f_highcut, fs, order, Is_detrend, Is_bandpass_filter); 

	Vel, Dis = acc_to_dis(acc[:, 0], acc_filtered_detrend); 

	Vel_original, Dis_original = acc_to_dis(acc[:, 0], acc[:, 1]);


	axial_label_font = FontProperties()
	axial_label_font.set_family('sans-serif')
	axial_label_font.set_style('normal')
	axial_label_font.set_weight('bold')
	# axial_label_font.set_size('x-large')
	axial_label_font.set_size(20)



	legend_label_font = FontProperties()
	legend_label_font.set_family('sans-serif')
	legend_label_font.set_style('normal')
	legend_label_font.set_weight('normal')
	# legend_label_font.set_size('large')
	legend_label_font.set_size(16)

	numbercol = 1; 

	fig = plt.figure(); 

	ax = fig.add_subplot(111)

	# ax.plot(acc[:, 0], acc[:, 1], '-k', label='original', linewidth= 2);

	# ax.plot(acc[:, 0], acc_filtered_detrend, ':r', label='filtered detrend', linewidth= 2);


	# ax.plot(f_original, fft_var_mag_original, '-r', label='original', linewidth= 3);

	# ax.plot(f_filtered, fft_var_mag_filtered, '-k', label='filtered', linewidth= 1.5);


	ax.plot(acc[:, 0], Dis_original, '-r', label='baseline corrected', linewidth= 2);

	ax.plot(acc[:, 0], Dis, '-k', label='original', linewidth= 2);

	# plt.gca().set_xlim([0,0.5]); 

	plt.gca().get_xaxis().set_tick_params(direction='in',labelsize='x-large')

	plt.gca().get_yaxis().set_tick_params(direction='in',labelsize='x-large')

	plt.xlabel('Time [s]', fontproperties=axial_label_font); 

	plt.ylabel('Dis. [$m$]', fontproperties=axial_label_font); 

	plt.grid(True); 

	plt.legend(ncol= numbercol, loc='lower left', prop=legend_label_font); 

	plt.savefig('baseline_correction.pdf', bbox_inches='tight'); 

	plt.show(); 






############################################################# 
############# Test region for Butterworth Bandpass ##########
############################################################# 

# fs = 5000.0

# lowcut = 500.0

# highcut = 1250.0

# T = 0.05

# nsamples = T * fs

# t = np.linspace(0, T, nsamples, endpoint=False)

# a = 0.02

# f0 = 600.0

# x = 0.1 * np.sin(2 * np.pi * 1.2 * np.sqrt(t))

# x += 0.01 * np.cos(2 * np.pi * 312 * t + 0.1)

# x += a * np.cos(2 * np.pi * f0 * t + .11)

# x += 0.03 * np.cos(2 * np.pi * 2000 * t)

# plt.figure()

# plt.clf()

# plt.plot(t, x, label='Noisy signal')

# y = butter_bandpass_filter(x, lowcut, highcut, fs, order=4)

# plt.plot(t, y, label='Filtered signal (%g Hz)' % f0)

# plt.xlabel('time (seconds)')

# plt.hlines([-a, a], 0, T, linestyles='--')

# plt.grid(True)

# plt.axis('tight')

# plt.legend(loc='upper left')

# plt.show()

