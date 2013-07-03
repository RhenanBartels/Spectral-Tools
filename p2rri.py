#!/usr/bin/env python
# _*_ coding: utf-8 _*_
"""
    Code for detection and visualization of RR-intervals that come from
    Polar acquisition board.

"""

import sys
from numpy import array, where, diff, cumsum
import matplotlib.pyplot as plt


def detpolar(x_is, thr, freqs):
    """
        Perform the detection of RR-intervals from the signal generated
        by a Polar board

    """
    x_is[x_is >= thr] = 1
    x_is[x_is < thr] = 0
    ct1 = x_is[0:-1]
    ct2 = x_is[1::]
    ctrl = ct2 - ct1
    pos = where(ctrl == 1)[0] + 1
    rri = diff(pos) / float(freqs)
    t_is = cumsum(rri)
    return t_is, rri

if __name__ == '__main__':

    if 'linux' in sys.platform:
        plt.switch_backend('qt4Agg')
    try:
        sys.argv[1]
    except IndexError:
        print "You must add a file. $python p2rri.py filename.txt"
        sys.exit()

    if 'txt' not in sys.argv[1]:
        raise Exception("File must be a text file")

    DT = []

    with open(sys.argv[1]) as f:
        for ln in f:
            try:
                DT.append(float(ln.strip()))
            except ValueError:
                continue

    TIME, RR, = detpolar(array(DT), 1, 1000)

    with open("RRi_" + sys.argv[1], 'w') as f:
        for ln in RR:
            f.write(str(ln) + "\n")

    plt.plot(TIME, RR)
    plt.title("Tachogram")
    plt.xlabel("TIME (s)")
    plt.ylabel("RRi (ms)")
    plt.axis("tight")
    plt.show()
