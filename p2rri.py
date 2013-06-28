#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys
from numpy import array, where, diff, cumsum
import matplotlib.pyplot as plt

def detPolar(x, thr, fs):
    x[x >= thr] = 1
    x[x < thr] = 0
    ct1 = x[0:-1]
    ct2 = x[1:len(x)]
    ct = ct2- ct1
    pos =where(ct == 1)[0] + 1
    rri = diff(pos) / 1000.0
    t = cumsum(rri)
    return t, rri

if __name__ == '__main__':

    if 'linux' in sys.platform:
            plt.switch_backend('qt4Agg')
    try:
        sys.argv[1]
    except IndexError:
        print "You must add a file. $ python p2rri.py filename.txt"

    if 'txt' not in sys.argv[1]:
        raise Exception("File must be a text file")

    dt = []

    with open(sys.argv[1]) as f:
        for ln in f:
            try:
                dt.append(float(ln.strip()))
            except ValueError:
                continue

    time, rri, = detPolar(array(dt), 1, 1000)
    plt.plot(time, rri)
    plt.title("Tachogram")
    plt.xlabel("Time (s)")
    plt.ylabel("RRi (ms)")
    plt.axis("tight")
    plt.show()
