import os
import numpy as np
import pylab as pl


def get_shot(type='g', mlim=25.5):
    files  = {'BX': 'BXDrop.dat', 'u': 'MalkanDrop.dat', 'g': 'gDrop.dat', 'r': 'rDrop.dat', 'qso': 'qso.dat'}
    file   = files[type]

    root   = os.environ['LBGCMB']
    fpath  = root + '/dropouts/schechter/dat/' + file        ## qso, BXDrop, MalkanDrop, gDrop, rDrop. 
    
    ms, ns = np.loadtxt(fpath, unpack=True)                  ## per sq. deg.  

    index  = (np.abs(ms - mlim) == np.abs(ms - mlim).min())

    mm     = ms[index][0]
    nn     = ns[index][0]
    
    return nn 

if __name__ == '__main__':
    print('\n\nWelcome to get_shot.\n\n')

    mlim = 25.5

    for type in ['BX', 'u', 'g', 'r', 'qso']:
      mm, nn = get_shot(type=type, mlim=mlim)
    
      print(mlim, nn)

    print('\n\nDone.\n\n')
