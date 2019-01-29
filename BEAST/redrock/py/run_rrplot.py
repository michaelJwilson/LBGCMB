##  Find where you installed rrplot
import os

rrplot = None

for p in os.environ['PATH'].split(':'):
    rrplot = os.path.join(p, 'rrplot')

    if os.path.exists(rrplot):
        break

if rrplot is None:
    print('ERROR: unable to find rrplot in your $PATH')

else:
    print('Using ' + rrplot)

#- Input files
specfile = 'spectro/redux/dc17a2/spectra-64/172/17242/spectra-64-17242.fits'
rrfile   = 'spectro/redux/dc17a2/spectra-64/172/17242/rr-64-17242.h5'

RRH5=os.environ['RRH5'] 
ZBEST=os.environ['ZBEST'] 

#- Now actually run it
%run $rrplot --specfile $specfile --rrfile $rrfile
