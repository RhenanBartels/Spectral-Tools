#-------------------------------------------------------------------------------
# Name:        pwelch
# Purpose:     Estimates the Modified Periodogram (Welch Peridogram) of a Signal
#
# Author:      Rhenan Bartels Ferreira
#
# Created:     25/03/2013
# Copyright:   (c) Rhenan 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import numpy as np


def pwelch(data, D, overlap, fs, side='onesided'):
    """If data x[0], x[1], ..., x[N-1] of N samplesis divided into P segments of D
    samples each, with a shift os S samples between adjacents segments (S<=D),
    then the maximum number of segments P is given by the integer part of:
    P = (N-D)/S + 1
    The  periodogram estivative is defined by:
    Pxx(p) = (1/(U*D*T))*abs(np.fft.fft(x[n]))**2
    where p is the range of segments O <= p <= P - 1
    over the frequency  -1/2*T <= f <= 1/2*T, where X(f) is the DFT of the pth
    segment. U is the discrete-time window energy and T is the sample interval.
    X(f) = T*x[n]*np.exp(-1j*2*np.pi*f*n*T)
    T = 1/fs
    U = T*np.sum(w[n]**2)"""

    start = 0
    end = D

    T = 1.0 / fs  # T is sample interval
    trun_size = len(data)/D
    data = data[0:trun_size*D]

    S = D - overlap
    P = (len(data) - D) / S + 1
    U = (1.0 / D) * sum(np.hanning(D)**2)  # Proakis Monolakis Pag 911
    #U = T*sum(np.hanning(D)**2)  #Digital Spectral Analysis Marple

    data_temp = []
    data_dc_han = []
    data_dc2_han = []
    PSD = []

    for i in xrange(P):
        data_temp = data[start:end]
        data_dc_han = (data_temp - np.mean(data_temp)) * np.hanning(D)
        data_dc2_han = data_dc_han - np.mean(data_dc_han)
        #Xf = abs(np.fft.fft(data_dc2_han/fs))**2 #Digital Spectral Analysis Marple
        Xf = (abs(np.fft.fft(data_dc2_han)) / D)**2  # Proakis Monolakis Pag 911
        PSD.append(Xf / (U))  # Digital Spectral Analysis Marple
        #PSD.append(Xf/(D*U)) # Proakis Monolakis Pag 911
        #PSD.append(Xf/(D*U)) # Notes on Digital Signal Processing: Note 31
        start += S
        end += S

    if side == 'onesided':
        Pxx = 2*np.sum(PSD, axis=0) / P
        Pxx = Pxx[0:D / 2.0 + 1]
        f = np.linspace(0, fs / 2, D / 2.0 + 1)
    else:
        Pxx = np.sum(PSD, axis=0) / P
        f = np.linspace(0, fs, len(Pxx))

    return [f, Pxx]
