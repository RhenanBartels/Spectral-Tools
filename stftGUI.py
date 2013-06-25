#!/usr/bin/env python
# _*_ coding: utf-8 _*_
import pygtk
pygtk.require('2.0')
import gtk
import matplotlib.pyplot as plt
import Tkinter
import tkFileDialog
from numpy import array, linspace, mean, loadtxt, fft, cumsum, arange, interp,\
        hanning

plt.switch_backend('qt4Agg')
plt.ion()

class Base:
    def destroy(self, widget, data=None):
        gtk.main_quit()


    def load(self, widget):
        master = Tkinter.Tk()
        master.withdraw()
        Path = tkFileDialog.askopenfilename(title="Choose the Signal",
                filetypes=[("Text File", "*.txt")])
        self.data = loadtxt(Path)


    def stft(self, widget):
        D = int(self.textbox1.get_text().strip())
        O = int(self.textbox2.get_text().strip())
        FS = int(self.textbox3.get_text().strip())
        t = cumsum(self.data) / 1000.0
        tx = arange(t[0], t[-1], 1.0 / FS)
        data = interp(tx, t, self.data)
        start = 0
        stop = D
        win = hanning
        S = D - O
        U = (1.0 / D) * sum(win(D))**2
        P = int(len(data) - D) / S + 1
        Sxx = []
        Time = []
        for k in xrange(P):
            data_temp = data[start:stop]
            data_dc = (data_temp - mean(data_temp)) * win(D)
            data_dc2 = data_dc - mean(data_dc)
            Xf = abs(fft.fft(data_dc2 / U))**2
            Sxx.append(Xf)
            Time.append(start)
            start += S
            stop += S
        Sxx = array(Sxx)[:, 0: FS / 2.0 + 1]
        print range(0, int(FS / 2.0) + 1)
        F = linspace(0, FS / 2.0, D / 2.0 + 1)
        Time = array(Time)
        #Time = linspace(0, start, Sxx.shape[0])
        plt.plot(tx, data)

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
        self.textbox1.set_text("512")
        self.textbox2.set_text("256")
        self.textbox3.set_text("4")
        fixed = gtk.Fixed()
        fixed.put(self.button1, 15, 30)
        fixed.put(self.button2, 75, 30)
        fixed.put(self.textbox1, 15, 60)
        fixed.put(self.textbox2, 15, 90)
        fixed.put(self.textbox3, 15, 120)
        self.window.add(fixed)
        self.window.show_all()
        self.window.connect("destroy", self.destroy)
    def main(self):
        gtk.main()

if __name__ == "__main__":
    base = Base()
    base.main()
