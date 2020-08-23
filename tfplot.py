# pylint: disable=invalid-name, missing-function-docstring, missing-module-docstring, too-many-locals, fixme

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
#   The sampling frequency is set by FS (default is 44.1 kHz) and the plotting color is
#   set by the COLOR variable (default is 'b').
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

import math
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft

# Ported from Octave fftfilt.m

# Constants
OCTBIN = 100.
FFTSIZE = 2**18
LOGFACT = 2**(1. / OCTBIN)


def fftfilt(b, x):

    # Use FFT with the smallest power of 2 which is >= length (x) + length (b) - 1
    # as number of points ...
    c_x = x.size
    c_b = b.size
    N = int(np.power(2, np.ceil(np.log(c_x + c_b) / np.log(2)), dtype=np.float32))
    y = ifft(fft(x, N) * fft(b, N))
    # Final cleanups: Both x and b are real; y should also be
    return np.real(y)


def tfplots(data, Fs=44100, color='b', fract=3):

    logn = np.floor(np.log(Fs / 2) / np.log(LOGFACT))
    # logarithmic scale from 1 Hz to Fs/2
    logscale = np.power(LOGFACT, np.r_[:logn])

    # creating a half hanning window
    wl = data.size
    hann = np.hanning(wl * 2)
    endwin = hann[wl:2 * wl]
    tf = fft(data * endwin, FFTSIZE)

    magn = np.abs(tf[:math.floor(FFTSIZE / 2)])

    # creating 100th octave resolution log. spaced data from the lin. spaced
    # FFT data
    logmagn = np.empty(int(logn))  # TODO: Validate conversion
    fstep = Fs / np.float64(FFTSIZE)

    for k in range(logscale.size):
        start = np.round(logscale[k] / np.sqrt(LOGFACT) / fstep)
        start = np.maximum(start, 1)
        start = np.minimum(start, FFTSIZE / 2)
        stop = np.round(logscale[k] * np.sqrt(LOGFACT) / fstep)
        stop = np.maximum(stop, 1)
        stop = np.minimum(stop, FFTSIZE / 2)
        start = int(start)
        stop = int(stop)

        # averaging the power
        logmagn[k] = np.sqrt(np.mean(np.power(magn[start - 1:stop], 2)))

    # creating hanning window
    # fractional octave smoothing
    hl = int(2 * np.round(OCTBIN / fract))
    hh = np.hanning(hl)

    L = logmagn.size
    logmagn[L - 1:L + hl] = 0

    # Smoothing the log. spaced data by convonvling with the hanning window
    tmp = fftfilt(hh, np.power(logmagn, 2))
    smoothmagn = np.sqrt(tmp[int(hl / 2):int(hl / 2) + L] / hh.sum(axis=0))

    # plotting
    plt.semilogx(logscale, 20 * np.log10(smoothmagn), color)


def tfplot(data, Fs=44100, color='b', avg='comp'):

    logn = np.floor(np.log(Fs / 2) / np.log(LOGFACT))
    # logarithmic scale from 1 Hz to Fs/2
    logscale = np.power(LOGFACT, np.r_[:logn])

    # creating a half hanning window
    wl = data.size
    hann = np.hanning(wl * 2)
    endwin = hann[wl:2 * wl]
    tf = fft(data * endwin, FFTSIZE)
    compamp = tf[:math.floor(FFTSIZE / 2)]

    logmagn = np.empty(math.floor(logn))
    fstep = Fs / np.float64(FFTSIZE)

    for k in range(logscale.size):

        # finding the start and end positions of the logaritmic bin
        start = np.round(logscale[k] / np.sqrt(LOGFACT) / fstep)
        start = np.maximum(start, 1)
        start = np.minimum(start, FFTSIZE / 2)
        stop = np.round(logscale[k] * np.sqrt(LOGFACT) / fstep) - 1
        stop = np.maximum(stop, start)
        stop = np.maximum(stop, 1)
        stop = np.minimum(stop, math.floor(FFTSIZE / 2))
        start = int(start)
        stop = int(stop)

        # averaging the complex transfer function
        if avg == 'comp':
            logmagn[k] = np.abs(np.mean(compamp[start - 1:stop]))
        elif avg == 'abs':
            logmagn[k] = np.mean(np.abs(compamp[start - 1:stop]))
        elif avg == 'power':
            logmagn[k] = np.sqrt(
                np.mean(np.abs(np.power(compamp[start - 1:stop], 2))))

    # plotting
    plt.semilogx(logscale, 20 * np.log10(logmagn), color)


def debug_log_plot(x, y):
    fig = plt.figure()
    plt.title("Digital filter frequency response")
    fig.add_subplot(111)
    plt.semilogx(x, y, 'b')
    plt.ylabel('Amplitude (power)', color='b')
    plt.xlabel('Frequency (rad/sample)')
    plt.grid()
    plt.legend()
    plt.show()
