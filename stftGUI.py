#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import pygtk
pygtk.require('2.0')
import gtk
import matplotlib.pyplot as plt
from numpy import array, linspace, mean, loadtxt, fft, cumsum, arange, interp,\
        hanning, zeros, concatenate, log10, ndarray

plt.switch_backend('qt4Agg')
plt.ion()

class Base:
    def destroy(self, widget, data=None):
        gtk.main_quit()

    def load(self, widget):
        #Stuff for create a open file dialog
        dialog = gtk.FileChooserDialog("Choose the Signal",
                None,
                gtk.FILE_CHOOSER_ACTION_OPEN,
                (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                    gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_OK)
        gtkFilter = gtk.FileFilter()
        gtkFilter.set_name("Text Files")
        gtkFilter.add_pattern("*.txt")
        dialog.add_filter(gtkFilter)
        response = dialog.run()
        if not response == gtk.RESPONSE_CANCEL:
            self.data = loadtxt(dialog.get_filename())
        dialog.destroy()
        self.button2.show()


    def stft(self, widget):
        #Short-time Fourier Transform
        D = int(self.textbox1.get_text().strip())  #Segment Size
        O = int(self.textbox2.get_text().strip())  #Overlap
        nfft = int(self.textbox3.get_text().strip())  #NFFT
        FS = int(self.textbox4.get_text().strip())  #Sampling Frequency
        t = cumsum(self.data) / 1000.0  #Time created from RRi series
        tx = arange(t[0], t[-1], 1.0 / FS)  #Interpolated Time
        data = interp(tx, t, self.data)  #Interpolated RRi
        #Stuff for Error Handling.
        if not isinstance(self.data, ndarray):
            messagedialog = gtk.MessageDialog(self.window, 0,
                    gtk.MESSAGE_ERROR, gtk.BUTTONS_OK,
                    "File must be a array_like. Read operation cancelled.")
            messsagedialog.run()
            messagedialog.destroy()
        elif D != int(D) or D <= 0 or O != int(O) or O <= 0:
            messagedialog = gtk.MessageDialog(self.window, 0,
                    gtk.MESSAGE_ERROR, gtk.BUTTONS_OK,
                    "Segment and Overlap Size Must be an Interger. "
                    "Read operation cancelled.")
            messagedialog.run()
            messagedialog.destroy()
        elif FS != int(FS) or FS < 0:
            messagedialog = gtk.MessageDialog(self.window, 0,
                    gtk.MESSAGE_ERROR, gtk.BUTTONS_OK,
                    "FS Must be a Positive Interger. "
                    "Read operation cancelled.")
            messagedialog.run()
            messagedialog.destroy()
        elif nfft != int(nfft) or nfft < 0:
            messagedialog = gtk.MessageDialog(self.window, 0,
                     gtk.MESSAGE_ERROR, gtk.BUTTONS_OK,
                    "NFFT Must be an Interger. Read operation cancelled.")
            messagedialog.run()
            messagedialog.destroy()
        elif O > D:
            messagedialog = gtk.MessageDialog(self.window, 0,
                    gtk.MESSAGE_ERROR, gtk.BUTTONS_OK,
                    "Overlap Must be Smaller Than Segment Size. "
                    "Read operation cancelled.")
            messagedialog.run()
            messagedialog.destroy()
        elif D > len(data):
            messagedialog = gtk.MessageDialog(self.window, 0,
                    gtk.MESSAGE_ERROR, gtk.BUTTONS_OK,
                    "Segment Size greater than Signal Length. "
                   " Read operation cancelled.")
            messagedialog.run()
            messagedialog.destroy()
        else:
            if nfft != 0:  #Test if the user wants NFFT
                zplen = nfft - D
                zp = zeros(zplen)
                U =(1.0 / nfft) * sum(hanning(nfft)**2)  #Window Energy
            else:
                zp = []
                U = (1.0 / D) * sum(hanning(D))**2  #Window Energy
            start = 0
            stop = D
            win = hanning  #TODO:More window options (hamming,blackman,kaiser etc)
            S = D - O  #Shift Size
            P = int(len(data) - D) / S + 1  #Number of Iterations
            Sxx = []
            Time = []
            for k in xrange(P):
                data_temp = concatenate((data[start:stop], zp), axis=0)  #Concatenates with zero-padding
                data_dc = (data_temp - mean(data_temp)) * hanning(len(data_temp))  #DC removing and window multiplication
                data_dc2 = data_dc - mean(data_dc)  #Residual DC from window multiplication removing
                Xf = 2 * abs(fft.fft(data_dc2 / len(data_temp)))**2  #PSD
                Sxx.append(Xf / U)  #Window Energy Compensation
                start += S
                stop += S
            Sxx = array(Sxx)[:, 0: D / 2.0 + 1]  #Get Rid of Half of Array
            Sxx = 10 * log10(Sxx)  #Perhaps better visualization
            Sxx[Sxx < 0] = 0
            F = linspace(0, FS / 2.0, D / 2.0 + 1)  #Frequency Axis
            Time = linspace(0, tx[-1], Sxx.shape[0])
            fig =  plt.figure()
            ax1 = fig.add_subplot(2, 1, 1)
            ax1.plot(tx, data)
            ax2 = fig.add_subplot(2, 1, 2)
            ax2.pcolor(Time, F, Sxx.T)  #TODO: Improve Plot Type
            ax2.set_ylim((0, 0.5))
            ax1.set_xlim((0, tx[-1]))
            ax2.set_xlim((0, Time[-1]))
            ax1.set_ylabel('RRi (ms)')
            ax2.set_ylabel('PSD (ms^2/Hz)')
            ax2.set_xlabel('Time (s)')

            #plt.imshow(Sxx, aspect='auto', interpolation='nearest')
    def __init__(self):
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_position(gtk.WIN_POS_CENTER)
        self.window.set_size_request(200, 200)
        self.button1 = gtk.Button("Load")
        self.button2 = gtk.Button("STFT")
        self.button1.connect("clicked", self.load)
        self.button2.connect("clicked", self.stft)
        self.textbox1 = gtk.Entry()
        self.textbox2 = gtk.Entry()
        self.textbox3 = gtk.Entry()
        self.textbox4 = gtk.Entry()
        self.textbox1.set_text("512")
        self.textbox2.set_text("256")
        self.textbox3.set_text("0")
        self.textbox4.set_text("4")
        fixed = gtk.Fixed()
        fixed.put(self.button1, 15, 30)
        fixed.put(self.button2, 75, 30)
        fixed.put(self.textbox1, 15, 60)
        fixed.put(self.textbox2, 15, 90)
        fixed.put(self.textbox3, 15, 120)
        fixed.put(self.textbox4, 15, 150)
        self.window.add(fixed)
        self.window.show_all()
        self.button2.hide()
        self.window.connect("destroy", self.destroy)
    def main(self):
        gtk.main()

if __name__ == "__main__":
    base = Base()
    base.main()
