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


def pwelch(data, D, overlap, fs, nfft=None, side='onesided'):
    """
    Parameters
    ----------
    data: array_like, 1-D
        Signal.
    D: interger
        Segment size.
    overlap: interger
        overlap of adjancent segments.
   fs: integer
       sampling frequency
   nfft: integer
       Size of Zero-Padding
   side: str or None, optional
       Must be 'onesided' or 'twosided'. This determinates the length of the
       outputs. If 'onesided' then 'f' and 'Pxx' goes from 0 to fs / 2 + 1,
       else, 'f' and 'Pxx' goes from 0 to fs -1.

   Returns
   ------
   f: array_like, 1-D
       frequency content.
   Pxx: array_like, 1_D
       Power Spectral Density
   ---------------------------------------------------------------------------

    If data x[0], x[1], ..., x[N - 1] of N samples is divided into P segments
    of D samples each, with a shift os S samples between adjacents segments
    (S<=D), then the maximum number of segments P is given by the integer part
    of: P = (N - D) / S + 1. The  periodogram estivative is defined by:

                  Pxx(p) = (1 / (U * D * T)) * abs(X(f))**2

    Where p is the range of segments O <= p <= P - 1 over the frequency
    -1/2*T <= f <= 1/2*T, where X(f) is the DFT of the pth segment.
    U is the discrete-time window energy and T is the sample interval.
    X(f) =  T * x[n] * np.exp(-1j * 2 * np.pi * f * n * T)
    T = 1 / fs
    U = T * np.sum(w[n]**2)

    """
    if not isinstance(data, np.ndarray):
        raise Exception("data must be a array_like")
    elif D != int(D):
        raise Exception("D must be a integer")
    elif overlap != int(overlap):
        raise Exception("overlap must be a integer")
    elif fs != int(fs):
        raise Exception("fs must be a integer")
    elif not isinstance(side, str):
        raise Exception("fs must be a string. 'one-sided' or 'twosided'")
    elif nfft is not None and nfft != int(nfft):
        raise ValueError("nfft must be None or an integer")
    elif nfft is not None and nfft <= D:
        raise Exception("nfft must be greater than D")
    elif overlap >= D:
        raise Exception("overlap must be smaller than D")
    if len(data) < D:
        raise Exception("D must be smaller than the length of data")


    start = 0
    end = D

    T = 1.0 / fs  # T is sample interval
    trun_size = len(data)/D
    data = data[0:trun_size*D]

    S = D - overlap
    P = int((len(data) - D) / S) + 1
    U = (1.0 / D) * sum(np.hanning(D)**2)  # Proakis Monolakis Pag 911
    #U = T * sum(np.hanning(D)**2)  #Digital Spectral Analysis Marple

    data_temp = []
    data_dc_han = []
    data_dc2_han = []
    PSD = []

    if nfft is not None:
        zplen = nfft - D
        zp = np.zeros(zplen)
        U = (1.0 / nfft) * sum(np.hanning(nfft)**2)
    else:
        zp = []

    for p in xrange(P):
        data_temp = np.concatenate((data[start:end], zp), axis=0)
        data_dc_han = (data_temp - np.mean(data_temp)) *\
                       np.hanning(len(data_temp))
        data_dc2_han = data_dc_han - np.mean(data_dc_han)
        #Xf = abs(np.fft.fft(data_dc2_han/fs))**2 #Digital Spectral Analysis Marple
        Xf = (abs(np.fft.fft(data_dc2_han)) / len(data_temp))**2  # Proakis Monolakis Pag 911
        PSD.append(Xf / U)  # Digital Spectral Analysis Marple
        #PSD.append(Xf/(D*U)) # Proakis Monolakis Pag 911
        #PSD.append(Xf/(D*U)) # Notes on Digital Signal Processing: Note 31
        start += S
        end += S

    if side == 'onesided':
        Pxx = 2 * np.sum(PSD, axis=0) / P
        Pxx = Pxx[0:D / 2.0 + 1]
        f = np.linspace(0, fs / 2, D / 2.0 + 1)
    else:
        Pxx = np.sum(PSD, axis=0) / P
        f = np.linspace(0, fs, len(Pxx))

    return f, Pxx
