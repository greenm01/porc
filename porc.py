#!/usr/bin/python -OO
#
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
#   More details about the parallel filter can be found in the papers
#
#   Balazs Bank, "Perceptually Motivated Audio Equalization Using Fixed-Pole Parallel
#   Second-Order Filters", IEEE Signal Processing Letters, 2008.
#   http://www.acoustics.hut.fi/go/spl08-parfilt
#
#   Balazs Bank, "Direct Design of Parallel Second-order Filters for
#   Instrument Body Modeling", International Computer Music Conference,
#   Copenhagen, Denmark, Aug. 2007.
#   http://www.acoustics.hut.fi/go/icmc07-parfilt
#
#   For Mixed-Phase Compensation, see:
#   "Mixed Time-Frequency approach for Mulitpoint Room Response Equalization," by
#   Alberto Carini, et al.

# Python libs
import sys
import textwrap
import wave
from contextlib import closing
import struct

# Scipy, Numpy, and matplotlibs
import numpy as np
import scipy as sp
import scipy.io as sio
import scipy.signal as sig
from scipy.fftpack import ifft, fft
from scipy.interpolate import pchip
from scipy.io import wavfile
from scipy.signal import convolve as conv
from scipy.stats import kurtosis
from scipy.stats import norm as Gaussian
from scipy import nanstd
import matplotlib.pyplot as plt

# PORC source files
from parfiltid import parfiltid
from tfplot import tfplot, tfplots, debug_log_plot
from freqpoles import freqpoles

# Ignore warnings
import warnings; warnings.filterwarnings('ignore')

# MiniDSP's OpenDRC box likes 6144 taps

def rceps(x): 
	y = sp.real(ifft(sp.log(sp.absolute(fft(x)))))
	n = len(x) 
	if (n%2) == 1:
		ym = np.hstack((y[0], 2*y[1:n/2], np.zeros(n/2-1)))
	else:
		ym = np.hstack((y[0], 2*y[1:n/2], y[n/2+1], np.zeros(n/2-1)))
	ym = sp.real(ifft(sp.exp(fft(ym)))) 
	return (y, ym)
    
def parfilt(Bm, Am, FIR, x):
    y = np.zeros(x.size)
    for k in range(Am.shape[1]):
        y += np.ravel(sig.lfilter(Bm[:,k], Am[:,k], x))
    y += np.ravel(sig.lfilter(np.hstack([FIR]), np.hstack([1]), x))
    return y
    
# Normalize signal
def norm(y): return y/np.fabs(y).max()

def dB2Mag(dB):
    return 10**((dB)/20.)
 
# The Median Absolute Deviation along given axis of an array
# From statsmodels lib
def mad(a, c=Gaussian.ppf(3/4.), axis=0):  # c \approx .6745
    a = np.asarray(a)
    return np.median((np.fabs(a))/c, axis=axis)
    
def roomcomp(impresp, filter, target, ntaps, mixed_phase, opformat, trim, nsthresh, noplot):

  print("Loading impulse response")
  
  # Read impulse response
  Fs, data = wavfile.read(impresp)
  data = norm(np.hstack(data))


  if trim:
    print("Removing leading silence")
    for spos,sval in enumerate(data):
        if abs(sval)>nsthresh:
            lzs=max(spos-1,0)
            ld =len(data)
            print('Impulse starts at position ', spos, '/', len(data))
            print('Trimming ', float(lzs)/float(Fs), ' seconds of silence')
            data=data[lzs:len(data)] #remove everything before sample at spos
            break
		  
  print("\nSample rate = ", Fs)
  
  print("\nGenerating correction filter")

  ###
  ## Logarithmic pole positioning
  ###

  fplog = np.hstack((sp.logspace(sp.log10(20.), sp.log10(200.), 14.), sp.logspace(sp.log10(250.), 
             sp.log10(20000.), 13.))) 
  plog = freqpoles(fplog, Fs)

  ###
  ## Preparing data
  ###

  # making the measured response minumum-phase
  cp, minresp = rceps(data)

  # Impulse response
  imp = np.zeros(len(data), dtype=np.float64)
  imp[0]=1.0

  # Target
  outf = []
  db = []

  if target == 'flat':
    
    # Make the target output a bandpass filter
    Bf, Af = sig.butter(4, 30/(Fs/2), 'high')
    outf = sig.lfilter(Bf, Af, imp) 
    
  else:
    
    # load target file
    t = np.loadtxt(target)
    frq = t[:,0]; pwr = t[:,1]
    
    # calculate the FIR filter via windowing method
    fir = sig.firwin2(501, frq, np.power(10, pwr/20.0), nyq = frq[-1])	
    # Minimum phase, zero padding	
    cp, outf = rceps(np.append(fir, np.zeros(len(minresp) - len(fir))))
      
  ###
  ## Filter design
  ###

  #Parallel filter design
  (Bm, Am, FIR) = parfiltid(minresp, outf, plog)

  # equalized loudspeaker response - filtering the 
  # measured transfer function by the parallel filter
  equalizedresp = parfilt(Bm, Am, FIR, data)

  # Equalizer impulse response - filtering a unit pulse
  equalizer = norm(parfilt(Bm, Am, FIR, imp))

  # Windowing with a half hanning window in time domain
  han = np.hanning(ntaps*2)[-ntaps:]
  equalizer = han * equalizer[:ntaps]

  ###
  ## Mixed-phase compensation
  ## Based on the paper "Mixed Time-Frequency approach for Multipoint
  ## Room Rosponse Equalization," by A. Carini et al.
  ## To use this feature, your Room Impulse Response should have all
  ## the leading zeros removed.
  ###
  if mixed_phase is True:
    
    # prototype function
    hp = norm(np.real(equalizedresp))

    # time integration of the human ear is ~24ms
    # See "Measuring the mixing time in auditoria," by Defrance & Polack
    hop_size = 0.024
    samples = hop_size * Fs

    bins = np.int(np.ceil(len(hp) / samples))

    tmix = 0

    # Kurtosis method
    for b in range(bins):
      start = np.int(b * samples)
      end = np.int((b+1) * samples)
      k = kurtosis(hp[start:end])
      if k <= 0:
        tmix = b * hop_size
        break

    # truncate the prototype function
    taps = np.int(tmix*Fs)

    print("\nmixing time(secs) = ", tmix, "; taps = ", taps)
    
    if taps > 0:
      # Time reverse the array
      h = hp[:taps][::-1]
      # create all pass filter
      phase = np.unwrap(np.angle(h))
      H = np.exp(1j*phase)
      # convert from db to linear
      mixed = np.power(10, np.real(H)/20.0)
      # create filter's impulse response
      mixed = np.real(ifft(mixed))
      
      # convolve and window to desired length
      equalizer = conv(equalizer, mixed)
      equalizer = han * equalizer[:ntaps]
      
      #data = han * data[:ntaps]
      #eqresp = np.real(conv(equalizer, data))
    else:
      print("zero taps; skipping mixed-phase computation")
  if opformat in ('wav', 'wav24'):      
  # Write data
    wavwrite_24(filter, Fs, norm(np.real(equalizer)))
    print('\nOutput format is wav24')
    print('Output filter length =', len(equalizer), 'taps')
    print('Output filter written to ' + filter)
	
    print("\nUse sox to convert output .wav to raw 32 bit IEEE floating point if necessary,")
    print("or to merge left and right channels into a stereo .wav")
    print("\nExample: sox leq48.wav -t f32 leq48.bin")
    print("         sox -M le148.wav req48.wav output.wav\n")

  elif opformat == 'wav32':
    wavwrite_32(filter, Fs, norm(np.real(equalizer)))
    print('\nOutput format is wav32')
    print('Output filter length =', len(equalizer), 'taps')
    print('Output filter written to ' + filter)
    print("\nUse sox to convert output .wav to raw 32 bit IEEE floating point if necessary,")
    print("or to merge left and right channels into a stereo .wav")
    print("\nExample: sox leq48.wav -t f32 leq48.bin")
    print("         sox -M le148.wav req48.wav output.wav\n")
  elif opformat == 'bin':
    # direct output to bin avoids float64->pcm16->float32 conversion by going direct 
    #float64->float32
    f = open(filter, 'wb')
    norm(np.real(equalizer)).astype('float32').tofile(f)
    f.close()
    print('\nOutput filter length =', len(equalizer), 'taps')
    print('Output filter written to ' + filter)
  else:
    print('Output format not recognized, no file generated.')


  ###
  ## Plots
  ###
  if not noplot:
    data *= 500
    # original loudspeaker-room response
    tfplot(data, Fs, avg = 'abs')
    # 1/3 Octave smoothed
    tfplots(data, Fs, 'r')

    #tfplot(mixed, Fs, 'r')

    # equalizer transfer function
    tfplot(0.75*equalizer, Fs, 'g')
    # indicating pole frequencies
    plt.vlines(fplog, -2, 2, color='k', linestyles='solid')

    # equalized loudspeaker-room response
    tfplot(equalizedresp*0.01, Fs, avg = 'abs')
    # 1/3 Octave smoothed
    tfplots(equalizedresp*0.01, Fs, 'r')

    # Add labels
    # May need to reposition these based on input data
    plt.text(325,30,'Unequalized loudspeaker-room response')
    plt.text(100,-15,'Equalizer transfer function')
    plt.text(100,-21,'(Black lines: pole locations)')
    plt.text(130,-70,'Equalized loudspeaker-room response')

    a = plt.gca()
    a.set_xlim([20, 20000])
    a.set_ylim([-80, 80])
    plt.ylabel('Amplitude (dB)', color='b')
    plt.xlabel('Frequency (Hz)')
    plt.grid()
    plt.legend()
    plt.show()

def wavwrite_24(fname, fs, data):
    data_as_bytes = (struct.pack('<i', int(samp*(2**23-1))) for samp in data)
    with closing(wave.open(fname, 'wb')) as wavwriter:
        wavwriter.setnchannels(1)
        wavwriter.setsampwidth(3)
        wavwriter.setframerate(fs)
        for data_bytes in data_as_bytes:
            wavwriter.writeframes(data_bytes[0:3])
            
def wavwrite_32(fname, fs, data):
    data_as_bytes = (struct.pack('<i', int(samp*(2**31-1))) for samp in data)
    with closing(wave.open(fname, 'wb')) as wavwriter:
        wavwriter.setnchannels(1)
        wavwriter.setsampwidth(4)
        wavwriter.setframerate(fs)
        for data_bytes in data_as_bytes:
            wavwriter.writeframes(data_bytes[0:4])
			
def main():
    
	print()

	mtxt = textwrap.dedent('''\
	Python Open Room Correction (PORC), version 0.1
	Copyright (c) 2012 Mason A. Green
	Based on the work of Dr. Balazs Bank
	''')

	bye = textwrap.dedent('''
	Example:
	./porc -t b&k.txt -n 8000 l48.wav leq48.bin
		
	See the README for detailed instructions
	''')

	import argparse
	from argparse import RawTextHelpFormatter

	parser = argparse.ArgumentParser(description = mtxt, epilog=bye, formatter_class=RawTextHelpFormatter)  

	# Positionals
	parser.add_argument('impresp', metavar='I', type=str, help='mesaured impulse response')
	parser.add_argument('filter', metavar='F', type=str, help='output filter file name')

	# Options
	parser.add_argument("-t", dest="target", default='flat',
					  help="target curve", metavar="FILE")
	parser.add_argument("-n", dest="ntaps", default = 6144,
					  help="filter length, in taps. Default = len(input)", type=int)
	parser.add_argument('--mixed', action='store_true', default = False,
					  help="Implement mixed-phase compensation. see README for details") 
	parser.add_argument("-o", dest="opformat", default = 'bin',
					help="Output file type, default bin optional wav", type=str)   
	parser.add_argument("-s", dest="nsthresh", default = 0.005,
					help="Normalized silence threshold. Default = 0.05", type=float)                    
	parser.add_argument('--trim', action='store_true', default = False,
					help="Trim leading silence")
	parser.add_argument('--noplot', action='store_true', default = False,
					help="Do not polt the filter") 	 					  

	args = parser.parse_args()

	roomcomp(args.impresp, args.filter, args.target, args.ntaps, args.mixed, args.opformat, args.trim, args.nsthresh, args.noplot)

if __name__=="__main__":
    main()  
