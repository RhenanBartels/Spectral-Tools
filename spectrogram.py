# _*_ coding: utf-8 _*_

from easygui import multenterbox
import Tkinter as Tk
import tkFileDialog
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import loadtxt, hanning, hamming, mean, linspace, array, interp,\
cumsum, arange
from numpy.fft import fft

def stft(data, segment, overlap, window, nfft, fs):
    shift = segment - overlap
    P = int(len(data) - segment) / shift + 1
    if window == "Hanning" or window == "hanning":
        win = hanning
    elif window == "Hamming" or window == "hamming":
        win = hamming
    else:
        raise Exception("Window not found")
    Sxx = []
    t_index = []
    U = sum(win(segment)**2) / segment
    start = 0
    stop = segment

    for k in xrange(P):
        data_temp = data[start:stop]
        data_dc = (data_temp - mean(data_temp)) * win(segment)
        data_dc2 = data_dc - mean(data_dc)
        xf = abs(fft(data_dc2) / segment )**2
        Sxx.append(xf / U)
        t_index.append(start)
        start += shift
        stop += shift

    Sxx = array(Sxx)
    Sxx = array([Sxx[i][0:segment / 2 + 1] for i in xrange(Sxx.shape[0])])
    t_index = array(t_index)
    F = linspace(0, fs / 2, segment / 2.0 + 1)
    return t_index, F, Sxx

if __name__ == '__main__':
    master = Tk.Tk()
    master.withdraw()
    path = tkFileDialog.askopenfilename(title = "Choose the data", filetypes=[
        ("text file", "*.txt")])
    eF = open(path)
    if not eF.name.endswith('txt'):
        raise Exception("File must be a text")
    else:
        try:
             data = loadtxt(path)
        except ValueError:
             print "Ops! Something went wrong"


    msg = "STFT Parameters"
    title = "Input Parameters"
    fieldNames = ["Segment size", "Overlap", "Window", "NFFT", "FS"]
    fieldValues = ["512", "256", "Hanning", "1024", "4"]
    fieldValues = multenterbox(msg, title, fieldNames, fieldValues)
    segment = int(fieldValues[0])
    overlap = int(fieldValues[1])
    window = fieldValues[2]
    nfft = int(fieldValues[3])
    fs = int(fieldValues[4])
    time = cumsum(data) / 1000.0
    tx = arange(time[0], time[-1], 1.0 / fs)
    data = interp(tx, time, data)
    t_index, F, Sxx = stft(data, segment, overlap, window, nfft, fs)
    plt.pcolor(t_index, F, Sxx.T)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
