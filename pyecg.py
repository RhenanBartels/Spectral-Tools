# _*_ coding: utf-8

from scipy.signal import butter, filtfilt
from easygui import multenterbox
import matplotlib.pyplot as plt
from numpy import diff, array, cumsum, where, arange, delete

def ecg2rri(x):
    global rri, t, ax2, pos_c
    pos_c = []
    plt.switch_backend('qt4Agg')
    plt.ion()
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
    t = cumsum(rri) / 1000.0
    t_ecg = arange(0, len(xf)) / Fs
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    ax1.plot(t_ecg, xf)
    ax1.plot(t_ecg[peaks], xf[peaks], 'g.-')
    ax2.plot(t, rri, 'k.-')
    fig.canvas.mpl_connect('button_press_event', onclick)
    return t, rri


def findMin(x, v):
    val = abs(v - x)
    pos = where(val == min(val))[0][0]
    print 'Position = %d, Value = %f' %(pos, v[pos])

    return val, pos


def onclick(event):
    global rri, t, ax2, poc_c
    xx = event.xdata
    val, pos = findMin(xx, t)
    pos_c.append(pos)
    if pos not in pos_c:
        rri = delete(rri, pos)
    t = cumsum(rri) / 1000.0
    ax2.plot(t, rri, 'k.-')
