#!/usr/bin/env python

"""
Extract PSF shape parameters vs. wavelength for use with quicksim.

Stephen Bailey, LBL
June 2014
"""

import sys
import os
import hashlib

import numpy as N
from scipy.interpolate import InterpolatedUnivariateSpline
import fitsio
import yaml

def calc_fwhm(x, y):
    """Return FWHM calculated by fitting a spline to y vs. x"""
    y = N.array(y)
    sp = InterpolatedUnivariateSpline(x, y-y.max()/2)
    lo, hi = sp.roots()
    return hi-lo

def img_fwhm(img):
    """Return FWHM(x), FWHM(y) of input 2D image"""
    x = N.arange(img.shape[1])
    y = N.arange(img.shape[0])
    fwhm_x = calc_fwhm(x, img.sum(axis=0))
    fwhm_y = calc_fwhm(y, img.sum(axis=1))

    return fwhm_x, fwhm_y

def calc_neff(img, pixsize):
    """Return effective number of cross-dispersion pixels for this psf spot"""
    #- Sum to x-axis and rebin to pixsize
    xpsf = img.sum(axis=0).reshape((pixsize, pixsize)).sum(axis=1)
    return N.sum(xpsf)**2 / N.sum(xpsf**2)

def quicksim_input_data(psffile, ww, ifiber=100):
    assert 0 <= ifiber < 500

    #- Read input data
    #- spots[i,j] is a 2D PSF spot sampled at
    #- slit position spotpos[i] and wavelength spotwave[j].
    #- Fiber k is located on the slit at fiberpos[k].
    spots = fitsio.read(psffile, 'SPOTS')
    spotwave = fitsio.read(psffile, 'SPOTWAVE')
    spotpos = fitsio.read(psffile, 'SPOTPOS')
    fiberpos = fitsio.read(psffile, 'FIBERPOS')
    hdr = fitsio.read_header(psffile)

    nwave = len(spotwave)
    npos = len(spotpos)
    nfiber = len(fiberpos)

    pixsize = int(round(float(hdr['CCDPIXSZ']) / hdr['CDELT1']))

    #- Measure the FWHM of the spots in x and y
    spot_fwhm_x = N.zeros((npos, nwave))
    spot_fwhm_y = N.zeros((npos, nwave))
    spot_neff = N.zeros((npos, nwave))
    for i in range(npos):
        for j in range(nwave):
            fx, fy = img_fwhm(spots[i,j])
            spot_fwhm_x[i,j] = fx
            spot_fwhm_y[i,j] = fy
            spot_neff[i,j] = calc_neff(spots[i,j], pixsize)

    #- For each spot wavelength, interpolate to the location of ifiber
    fiber_fwhm_x = N.zeros(nwave)
    fiber_fwhm_y = N.zeros(nwave)
    fiber_neff = N.zeros(nwave)

    for j in range(nwave):
        spx = InterpolatedUnivariateSpline(spotpos, spot_fwhm_x[:, j])
        fiber_fwhm_x[j] = spx(fiberpos[ifiber])
        spy = InterpolatedUnivariateSpline(spotpos, spot_fwhm_y[:, j])
        fiber_fwhm_y[j] = spy(fiberpos[ifiber])

        spn = InterpolatedUnivariateSpline(spotpos, spot_neff[:, j])
        fiber_neff[j] = spn(fiberpos[ifiber])


    #- Interpolate onto ww wavelength grid
    spx = InterpolatedUnivariateSpline(spotwave, fiber_fwhm_x)
    fwhm_x = spx(ww)
    spy = InterpolatedUnivariateSpline(spotwave, fiber_fwhm_y)
    fwhm_y = spy(ww)

    #- Convert fwhm units from spot pixels to CCD pixels
    #- Use units propagated from original spots calculations, not desi.yaml
    fwhm_x /= pixsize
    fwhm_y /= pixsize

    #- Final Neff sampled on same wavelength grid
    spn = InterpolatedUnivariateSpline(spotwave, fiber_neff)
    neff = spn(ww)

    #- Angstroms per row
    from numpy.polynomial.legendre import Legendre
    ycoeff, yhdr = fitsio.read(psffile, 'YCOEFF', header=True)
    domain = (yhdr['WAVEMIN'], yhdr['WAVEMAX'])
    y = Legendre(ycoeff[ifiber], domain=domain)(ww)
    ang_per_row = N.gradient(ww) / N.gradient(y)

    #- Convert fwhm_y from pixels to Angstroms
    fwhm_y *= ang_per_row

    data = N.rec.fromarrays([ww, fwhm_y, fwhm_x, neff, ang_per_row],
        names="wavelength,fwhm_wave,fwhm_spatial,neff_spatial,angstroms_per_row")

    return data

#-------------------------------------------------------------------------
desi = yaml.load(open(os.getenv('DESIMODEL')+'/data/desi.yaml'))

import argparse

parser = argparse.ArgumentParser(usage = "%prog [options]")
parser.add_argument("-o", "--output", action='store',  help="output fits file")
opts = parser.parse_args()

if opts.output is None:
    opts.output = os.path.join(os.getenv('DESIMODEL'), 'data', 'specpsf', 'psf-quicksim.fits')

clobber = True
for camera in ('b', 'r', 'z'):
    psffile = '{}/data/specpsf/psf-{}.fits'.format(os.getenv('DESIMODEL'), camera)
    psfsha1 = hashlib.sha1(open(psffile, mode='rb').read()).hexdigest()
    psfhdr = fitsio.read_header(psffile)

    wavemin = psfhdr['WMIN_ALL']
    wavemax = psfhdr['WMAX_ALL']

    #- The final FWHM grid is interpolated on 0.5 Angstrom grid
    dw = 0.5
    ww = N.arange(wavemin, wavemax+dw/2, dw)

    data = quicksim_input_data(psffile, ww)

    #- output header
    hdr = list()
    ### hdr.append(dict(name='EXTNAME', value='QUICKSIM-'+camera.upper(), comment='QuickSim params for camera '+camera))
    hdr.append(dict(name='PSFFILE', value=os.path.basename(psffile), comment='Input PSF file'))
    hdr.append(dict(name='PSFSHA1', value=psfsha1, comment='SHA1 checksum input PSF'))
    hdr.append(dict(name='WMIN_ALL', value=wavemin, comment='Starting wavelength [Angstroms]'))
    hdr.append(dict(name='WMAX_ALL', value=wavemax, comment='Last wavelength [Angstroms]'))
    hdr.append(dict(name='WAVEUNIT', value='Angstrom', comment='Wavelengths in Angstroms'))

    hdr.append(dict(name='TUNIT1', value='Angstrom', comment='Wavelength'))
    hdr.append(dict(name='TUNIT2', value='Angstrom', comment='Wavelength dispersion FWHM [Angstrom]'))
    hdr.append(dict(name='TUNIT3', value='pixel', comment='Cross dispersion FWHM [pixel]'))
    hdr.append(dict(name='TUNIT4', value='pixel', comment='Effective number of cross-dispersion pixels'))
    hdr.append(dict(name='TUNIT5', value='Angstrom/pixel', comment='Angstroms per CCD pixel row'))

    extname = 'QUICKSIM-'+camera.upper()
    fitsio.write(opts.output, data, header=hdr, clobber=clobber, extname=extname)
    clobber = False

#--- DEBUG ---
# import IPython
# IPython.embed()
#--- DEBUG ---
