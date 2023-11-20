# Python Open Room Correction (PORC)
# Copyright (c) 2012 Mason A. Green
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# TFPLOT - Smoothed transfer fucntion plotting
#   TFPLOTS(IMPRESP,COLOR, Fs, FRACT)
#   Logarithmic transfer function plot from impluse response IMPRESP. 
#   A half hanning window is applied before a 2^18 point FFT, then the data is colleced
#   into logaritmically spaced bins and the average power is computed for
#   each bin (100/octave). Then this is power-smoothed by a hanning window, where
#   FRACT defines the fractional-octave smoothing (default is 3, meaning third-octave).
#   The length of the smoothing hanning window is the double compared to the distance
#   defined by FRACT.
#   The sampling frequency is set by FS (default is 44.1 kHz) and the plotting color is set by the COLOR variable
#   (default is 'b').
#
#   TFPLOT(IMPRESP, COLOR, FS, OCTBIN, AVG)
#   Logarithmic transfer function plot from impluse response IMPRESP. 
#   A half hanning window is applied before a 2^18 point FFT, then the 
#   data is colleced into logaritmically spaced bins and the average 
#   response is computed for each bin. OCTBIN sets the number of bins 
#   in one octave, the default is 100 (lower numbers mean more smoothing). 
#   The sampling frequency is set by FS (default is 44.1 kHz) and the 
#   plotting color is set by the COLOR variable (default is 'b').
#
#   If the AVG variable is set to 'power' then the power is averaged
#   in the logaritmic bin, if it is 'abs' then the absolute value. If the
#   AVG parameter is set to 'comp' or omitted, it averages the complex
#   magnitude (i.e., this is the default).
#
#   C. Balazs Bank, 2006-2007.

import numpy as np
import scipy as sp
import scipy.signal as sig
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft

# Ported from Octave fftfilt.m
def fftfilt(b, x):

	## Use FFT with the smallest power of 2 which is >= length (x) +
	## length (b) - 1 as number of points ...
	c_x = x.size
	c_b = b.size
	N = np.power(2, np.ceil(np.log(c_x+c_b)/np.log(2)), dtype=np.float32)
	y = ifft(fft(x, N)*fft(b, N))
	## Final cleanups: Both x and b are real; y should also be
	return np.real(y)
	
def tfplots(data, Fs = 44100, color = 'b', fract=3):

	octbin = 100.
	FFTSIZE = 2**18

	logfact = 2**(1./octbin)
	LOGN = np.floor(np.log(Fs/2)/np.log(logfact))
	# logarithmic scale from 1 Hz to Fs/2
	logscale = np.power(logfact, np.r_[:LOGN]) 

	# creating a half hanning window
	WL = data.size
	hann = sp.hanning(WL*2)
	endwin = hann[WL:2*WL]
	tf = fft(data*endwin, FFTSIZE)

	magn = np.abs(tf[:FFTSIZE/2])
	compamp = tf[:FFTSIZE/2]

	# creating 100th octave resolution log. spaced data from the lin. spaced FFT data
	logmagn = np.empty(LOGN)
	fstep = Fs/np.float64(FFTSIZE)
	
	for k in range(logscale.size):
		start = np.round(logscale[k]/np.sqrt(logfact)/fstep)
		start = np.maximum(start,1)
		start = np.minimum(start, FFTSIZE/2)
		stop = np.round(logscale[k]*np.sqrt(logfact)/fstep)
		stop = np.maximum(stop,1)
		stop = np.minimum(stop, FFTSIZE/2)
		# averaging the power
		logmagn[k] = np.sqrt(np.mean(np.power(magn[start-1:stop],2))) 

	# creating hanning window
	# fractional octave smoothing
	HL = 2 * np.round(octbin/fract)
	hh = sp.hanning(HL)
	
	L = logmagn.size
	logmagn[L-1:L+HL] = 0

	# Smoothing the log. spaced data by convonvling with the hanning window
	tmp = fftfilt(hh, np.power(logmagn,2))
	smoothmagn = np.sqrt(tmp[HL/2:HL/2+L]/hh.sum(axis=0))

	# plotting
	plt.semilogx(logscale, 20*np.log10(smoothmagn), color)

def tfplot(data, Fs = 44100, color = 'b', octbin = 100, avg = 'comp'):

	FFTSIZE=2**18

	logfact = 2**(1./octbin)
	LOGN = np.floor(np.log(Fs/2)/np.log(logfact))
	# logarithmic scale from 1 Hz to Fs/2
	logscale = np.power(logfact, np.r_[:LOGN]) 

	# creating a half hanning window
	WL = data.size
	hann = sp.hanning(WL*2)
	endwin = hann[WL:2*WL]
	tf = fft(data*endwin, FFTSIZE)
	compamp = tf[:FFTSIZE/2]

	logmagn = np.empty(LOGN)
	fstep = Fs/np.float64(FFTSIZE)
	
	for k in range(logscale.size):

		#finding the start and end positions of the logaritmic bin
		start = np.round(logscale[k]/np.sqrt(logfact)/fstep)
		start = np.maximum(start, 1);
		start = np.minimum(start, FFTSIZE/2)
		stop = np.round(logscale[k]*np.sqrt(logfact)/fstep)-1
		stop = np.maximum(stop, start)
		stop = np.maximum(stop, 1)
		stop = np.minimum(stop, FFTSIZE/2)

		#averaging the complex transfer function
		if avg == 'comp':
			logmagn[k] = np.abs(np.mean(compamp[start-1:stop]))
		elif avg == 'abs':
			logmagn[k] = np.mean(np.abs(compamp[start-1:stop]))
		elif avg == 'power':
			logmagn[k] = np.sqrt(np.mean(np.abs(np.power(compamp[start-1:stop],2))))

	# plotting
	plt.semilogx(logscale, 20*np.log10(logmagn), color)
	
def debug_log_plot(x, y):
	fig = plt.figure()
	plt.title("Digital filter frequency response")
	ax = fig.add_subplot(111)
	plt.semilogx(x, y, 'b')
	plt.ylabel('Amplitude (power)', color='b')
	plt.xlabel('Frequency (rad/sample)')
	plt.grid()
	plt.legend()
	plt.show()
