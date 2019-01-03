#!/usr/bin/env python

"""
Make plots to compare two different versions of desimodel

Stephen Bailey, LBL
July 2014
"""

import os, sys
import numpy as np
import pylab as P
import matplotlib.pyplot as plt
import fitsio

camcolor = dict(b='b', r='r', z='k')

def compare_throughput(dir1, dir2):
    P.figure()
    p0 = plt.subplot2grid((3,1), (0,0), rowspan=2)
    p1 = plt.subplot2grid((3,1), (2,0))
    for x in ('b', 'r', 'z'):
        d1 = fitsio.read(dir1+'/data/throughput/thru-'+x+'.fits')
        d2 = fitsio.read(dir2+'/data/throughput/thru-'+x+'.fits')
        
        w1 = d1['wavelength']
        w2 = d2['wavelength']
        
        t1 = d1['throughput']
        t2 = d2['throughput']
    
        p0.plot(w1, t1, '-', color=camcolor[x])
        p0.plot(w2, t2, '--', color=camcolor[x])
        
        p1.plot(w1, (t1-np.interp(w1, w2, t2))/t1, '-', color=camcolor[x])

    p0.set_xlim(3500, 10000)
    p0.set_ylim(0.0, 0.5)
    p0.set_ylabel('Throughput')
    p0.grid()

    p1.set_xlim(3500, 10000)
    ### p1.set_ylim(-0.5, 0.5)
    p1.set_xlabel('Wavelength [Angstroms]')
    p1.set_ylabel('Relative difference')
    p1.grid()
    
def compare_fiberloss(dir1, dir2):
    pass

#-------------------------------------------------------------------------
dir1, dir2 = sys.argv[1:3]

compare_throughput(dir1, dir2)

plt.show()
