#!/usr/bin/python
# -*- coding: utf-8 -*-
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
# PARFILTID - System identification in the form of second-order parallel filters for a given
# pole set.
#   [Bm,Am,FIRcoeff]=parfiltid(INPUT,OUTPUT,P,NFIR); identifies the second-order sections
#   [Bm,Am] and the coefficients of the FIR part (FIRcoeff) for a given
#   pole set P. The parameters are set such that, when the INPUT signal is
#   filtered by the parallel filter, it gives an output which is the closest
#   to the OUTPUT vector in the LS sense. The number of taps in the parallel
#   FIR filter is set by NFIR. The default is NFIR=1, in this case FIRcoeff
#   is a simple gain. The only difference from the PARFILTDES function is that
#   now the input can be arbitrary, and not just a unit pulse as for filter design.
#
#   The Bm and Am matrices are containing the [b0 b1]' and [1 a0 a1]'
#   coefficients for the different sections in their columns. For example,
#   Bm(:,3) gives the [b0 b1]' parameters of the third second-order
#   section. These can be used by the filter command separatelly (e.g., by
#   y=filter(Bm(:,3),Am(:,3),x), or by the PARFILT command.
#
#   Note that this function does not support pole multiplicity, so P should
#   contain each pole only once.
#
#   More details about the parallel filter can be found in the papers
#
# ....Balazs Bank, "Perceptually Motivated Audio Equalization Using Fixed-Pole Parallel
#   Second-Order Filters", IEEE Signal Processing Letters, 2008.
#   http://www.acoustics.hut.fi/go/spl08-parfilt
#
#   Balazs Bank, "Direct Design of Parallel Second-order Filters for
#   Instrument Body Modeling", International Computer Music Conference,
#   Copenhagen, Denmark, Aug. 2007.
#   http://www.acoustics.hut.fi/go/icmc07-parfilt
#
#   C. Balazs Bank, Helsinki University of Technology, 2007.

# pylint: disable=invalid-name

import numpy as np
import scipy as sp
import scipy.signal as sig


def parfiltid(
    input,
    out,
    p,
    NFIR=1,
):

    # We don't want to have any poles in the origin; For that we have the parallel FIR part.
    # Remove nonzeros

    p = p[p.nonzero()]

    # making the filter stable by flipping the poles into the unit circle

    for k in range(p.size):
        if abs(p[k]) > 1:
            p[k] = 1.0 / np.conj(p[k])

    # Order it to complex pole pairs + real ones afterwards

    p = np.sort_complex(p)

    # in order to have second-order sections only (i.e., no first order)

    pnum = len(p)  # number of poles
    ppnum = 2 * np.floor(pnum / 2)  # the even part of pnum
    ODD = 0

    # if pnum is odd

    if pnum > ppnum:
        ODD = 1

    OUTL = len(out)
    INL = len(input)

    # making input the same length as the output

    if INL > OUTL:
        input = input[:OUTL]

    if INL < OUTL:
        input = np.hstack([input, np.zeros(OUTL - INL,
                                           dtype=np.float64)])

    L = OUTL

    # Allocate memory

    M = np.zeros((input.size, p.size + NFIR), dtype=np.float64)

    # constructing the modeling signal matrix....

    for k in range(0, int(ppnum), 2):  # second-order sections

        # impluse response of the two-pole filter

        resp = sig.lfilter(np.array([1]), np.poly(p[k:k + 2]), input)
        M[:, k] = resp

        # the response delayed by one sample

        M[:, k + 1] = np.hstack((0., resp[:L - 1]))

    # if the number of poles is odd, we have a first-order section

    if ODD:
        resp = sig.lfilter(np.array([1]), np.poly(p[-1]), input)
        M[:, pnum - 1] = resp

    # parallel FIR part

    for k in range(0, NFIR):
        M[:, pnum + k] = np.hstack([np.zeros(k, dtype=np.float64),
                                    input[:L - k + 1]])

    y = out

    # Looking for min(||y-M*par||) as a function of par:
    # least squares solution by equation solving

    mconj = M.conj().T
    A = np.dot(mconj, M)
    b = np.dot(mconj, y)
    par = np.linalg.solve(A, b)

    # print (np.dot(A, par) == b).all()

    # Allocate memory

    size = int(np.ceil(ppnum / 2))
    Am = np.zeros((3, size), dtype=np.float64)
    Bm = np.zeros((2, size), dtype=np.float64)

    # constructing the Bm and Am matrices

    for k in range(0, size):
        Am[:, k] = np.poly(p[2 * k:2 * k + 2])
        Bm[:, k] = np.hstack(par[2 * k:2 * k + 2])

    # we extend the first-order section to a second-order one by adding zero
    # coefficients

    if ODD:
        Am = np.append(Am, np.vstack(np.hstack([np.poly(p[pnum]),
                                                0.])), 1)
        Bm = np.append(Bm, np.vstack([par[pnum], 0.]), 1)

    FIR = []

    # constructing the FIR part

    if NFIR > 0:
        FIR = np.hstack(par[pnum:pnum + NFIR])

    return (Bm, Am, FIR)
