#!/usr/bin/python
# -*- coding: utf-8 -*-
# pylint: disable=invalid-name, missing-module-docstring, missing-function-docstring

# FREQPOLES - Pole set generation for the PARFILTDES function.
#   [P]=freqpoles(FR,Fs) creates poles from the given frequency vector
#   FR with sampling frequency Fs such that the frequency responses
#   of the second-order basis functions cross approximatelly at their
#   -3 dB point.
#   (See the equations in my 128th AES Convention paper.)
#
#   [P]=freqpoles(FR,Fs,Q) creates poles from the given frequency vector
#   FR with sampling frequency Fs but it uses the quality factors
#   given in Q vector to set the radius of the poles. (-3dB bandwidth
#   is DeltaF=FR/Q.)
#
#   C. Balazs Bank, 2010.

import numpy as np


def freqpoles(fr, Fs=44100, Q=None):

    # discrete pole frequencies

    wp = 2 * np.pi * fr / Fs

    dwp = np.empty(wp.shape)
    if Q is None:
        pnum = wp.size
        dwp = np.zeros(pnum)
        for k in range(1, pnum - 1):
            dwp[k] = (wp[k + 1] - wp[k - 1]) / 2
        dwp[0] = wp[1] - wp[0]
        dwp[pnum - 1] = wp[pnum - 1] - wp[pnum - 2]
    else:
        dwp = wp / Q

    p = np.exp(-dwp / 2) * np.exp(1j * wp)

    return np.hstack([p, p.conj()])
