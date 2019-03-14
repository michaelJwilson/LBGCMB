import os
import numpy as np
import pylab as pl


def get_shot(type='gDrop', mlim=25.5):
    root   = os.environ['LBGCMB']
    fpath  = root + '/dropouts/schechter/dat/%s.dat' % type  ## qso, BXDrop, MalkanDrop, gDrop, rDrop. 
    
    ms, ns = np.loadtxt(fpath, unpack=True)                  ## per sq. deg.  

    index  = (np.abs(ms - mlim) == np.abs(ms - mlim).min())

    mm     = ms[index][0]
    nn     = ns[index][0]
    
    return mm, nn 

if __name__ == '__main__':
    print('\n\nWelcome to get_shot.\n\n')

    mlim = 25.5

    for type in ['BXDrop', 'MalkanDrop', 'gDrop', 'rDrop']:
      mm, nn = get_shot(type=type)
    
      print(mm, nn)

    print('\n\nDone.\n\n')
