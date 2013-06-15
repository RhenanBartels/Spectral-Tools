# _*_ coding: utf-8

from scipy.signal import butter, filtfilt
from easygui import multenterbox
import matplotlib.pyplot as plt
from numpy import diff, array, cumsum, where

def ecg2rri(x):
    msg = "RRi Detection Parameters"
    title = "Parameters Dialog"
    fieldNames = ["Threshold", "Refractory Period", "Low Cuttof Freq.",
                  "Upper Cuttof Freq.", "Sampling Frequency"]
    fieldValues = ["0.5", "200", "5", "40", "1000"]
    fieldValues = multenterbox(msg, title, fieldNames, fieldValues)
    fieldValues = [float(k) for k in fieldValues]
    thr = fieldValues[0]
    Fs = fieldValues[4]
    lC = fieldValues[2] / (0.5 * Fs)  #normalized lower cutoff frequency.
    uC = fieldValues[3] / (0.5 * Fs)  #normalized upper cutoff frequency.
    B, A = butter(4, [lC, uC], 'band')
    xf = filtfilt(B, A, x)  #filtered ecg.
    xd = diff(xf)  #first derivative of the ecg.
    peaks = array([peaks + 1 for peaks in xrange(len(xd)) if xd[peaks] > 0 and
        xd[peaks + 1] < 0 and xf[peaks] > thr or xd[peaks] == 0])  #find RR peaks above threshold.

    rri = diff(peaks)  # RRi in miliseconds
    t = cumsum(rri)
    return t, rri, peaks

def findMin(x, y):
    val = []
    pos = []
    val = [min(abs(tt - x)) for tt in y]
    pos = [where(y == comp)[0][0] for comp in val]
    return val, pos


    return val, pos

def renderPlot(ax1, ax2, x1, x2):
    ax1.plot(x1)
    ax2.plot(x2)

def onclick(event):
    global rr_ind
    rr_ind.append(event.x)

def manualEdit(t_ecg, ecg, t_rri, rri):
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2,1, 2)
    ax1.plot(t_ecg, ecg)
    ax2.plot(t_rri, rri)
    rr_ind = []
    ev = fig.canvas.mpl_connect('button_press_event', onclick)
    vp, mp = findMin(t_rri, t_ecg[rr_ind])

