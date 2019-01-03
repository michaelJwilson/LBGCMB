#!/usr/bin/env python
# -*- coding: utf-8 -*-
from sys import exit
#
#
#
def main():
    """
    Convert simulated DESI spectrograph PSF spots into Specter PSF format.

    Spots and their CCD (x,y) location are provided on a grid of slit positions
    and wavelengths.  Fiber number and CCD x position increase with slit position;
    CCD y position increases with wavelength.  These spots and locations must
    be interpolated to the actual fiber positions on the slit and
    to arbitrary wavelengths.

    This code writes a Specter SpotGridPSF format to encode this information.

    Stephen Bailey, LBL
    September 2013
    """

    import sys
    import os
    import numpy as N
    from scipy import ndimage  #- for center_of_mass and shift
    from numpy.polynomial.legendre import Legendre
    import fitsio
    import yaml

    #- Load options
    import argparse
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    # parser.add_argument("-p", "--prefix", action='store',  help="input psf files prefix, including path")
    parser.add_argument("-o", "--outpsf", action='store',  help="output PSF file",
                        default='psf-blat.fits')
    # parser.add_argument("-t", "--throughput", action='store',  help="input throughput file to embed with PSF")
    parser.add_argument("-d", "--debug", action="store_true",  help="start ipython prompt when done")
    parser.add_argument("-c", "--camera", action='store',  help="camera: b, r, or z")
    parser.add_argument('spotfiles', action='store', help='Input spot files', narg='+')

    opts = parser.parse_args()

    if len(opts.spotfiles) == 0:
        print("ERROR: no input spot files given", file=sys.stderr)
        return 1

    #- Read DESI parameters
    params = yaml.load(open(os.getenv('DESIMODEL')+'/data/desi.yaml'))

    #- Get dimensions from first spot file
    hdr = fitsio.read_header(spotfiles[0])
    SpotPixelSize = hdr['PIXSIZE']  #- PSF spot pixel size in mm

    #- Hardcode spectrograph and CCD dimensions
    # CcdPixelSize = 0.015   #- CCD pixel size in mm
    # FiberSpacing = 0.230   #- center-to-center spacing in mm
    # GroupSpacing = 0.556   #- center-to-center group gap in mm
    # FibersPerGroup = 25
    # GroupsPerCcd = 20
    # NumFibers = 500
    # NumPixX = 4096
    # NumPixY = 4096
    # nspec = FibersPerGroup * GroupsPerCcd

    #- CCD pixel size in mm
    CcdPixelSize = params['ccd'][opts.camera]['pixsize'] / 1000.0  #- um -> mm

    #- center-to-center fiber spacing in mm on slit
    FiberSpacing = params['spectro']['fiber_spacing']

    #- center-to-center fiber group gap in mm on slit
    GroupSpacing = params['spectro']['fiber_group_spacing']

    FibersPerGroup = params['spectro']['fibers_per_group']
    GroupsPerCcd = params['spectro']['groups_per_ccd']
    NumFibers = params['spectro']['nfibers']
    NumPixX = params['ccd'][opts.camera]['npix_x']
    NumPixY = params['ccd'][opts.camera]['npix_y']
    nspec = FibersPerGroup * GroupsPerCcd

    #- Determine grid of wavelengths and fiber positions for the spots
    #- Use set() to get unique values, then convert to sorted array
    #- spotgrid maps (fiberpos, wavelength) -> filename
    print("Determining wavelength and slit position grid")
    wavelength = set()
    spotpos = set()
    spotgrid = dict()
    for filename in spotfiles:
        hdr = fitsio.read_header(filename)
        w = hdr['WAVE']*10      #- Wavelength [nm -> AA]
        p = hdr['FIBER']        #- Fiber slit position [mm]
        p = -p      #- Swap slit axis orientation to match CCD x
        wavelength.add(w)       #- Wavelength nm -> AA
        spotpos.add(p)
        spotgrid[(p,w)] = filename

    #- Wavelengths and slit positions of spots in grid
    wavelength = N.array( sorted(wavelength) )
    spotpos = N.array( sorted(spotpos) )

    #- Load grid of spots, and the x,y CCD pixel location of those spots
    print("Reading spots")
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    np = len(spotpos)
    nw = len(wavelength)
    spots = N.zeros( (np, nw, ny, nx), dtype=N.float32 )
    spotx = N.zeros( (np, nw), dtype=N.float32 )
    spoty = N.zeros( (np, nw), dtype=N.float32 )
    for i, p in enumerate(spotpos):
        for j, w in enumerate(wavelength):
            pix = fitsio.read(spotgrid[(p,w)])
            hdr = fitsio.read_header(spotgrid[(p,w)])

            #- Shift spot to center of image
            #- NOTE: uses spline interpolation, not sinc interpolation
            npy, npx = pix.shape
            yc,xc = ndimage.center_of_mass(pix)
            xmid = (pix.shape[1]-1)/2.0
            ymid = (pix.shape[0]-1)/2.0
            dx = xmid - xc
            dy = ymid - yc
            spots[i,j] = ndimage.shift(pix, (dy,dx))

            #- Reference pixel in FITS file
            xref = hdr['CRPIX1']-1
            yref = hdr['CRPIX2']-1

            #- Location of centroid on CCD in mm from center
            spotx[i,j] = hdr['CRVAL1'] + (xmid-xref+dx)*hdr['CDELT1']
            spoty[i,j] = hdr['CRVAL2'] + (ymid-yref+dy)*hdr['CDELT2']

    #- Convert spotx, spoty to pixel units instead of mm
    spotx = spotx/CcdPixelSize + NumPixX/2
    spoty = spoty/CcdPixelSize + NumPixY/2

    #- Map location of each fiber along the slit
    ifiber = N.arange(NumFibers).astype(int)
    ngaps = ifiber / FibersPerGroup    #- Number of gaps prior to fiber ifiber
    fiberpos = ifiber*FiberSpacing + ngaps*(GroupSpacing - FiberSpacing)
    fiberpos -= N.mean(fiberpos)

    #-----
    #- Determine range of wavelengths to fit
    #- Fit Legendre polynomials and extrapolate to CCD edges
    wmin = wavelength[0]
    wmax = wavelength[-1]
    for i in range(np):
        poly = Legendre.fit(spoty[i], wavelength, deg=5, domain=(0, NumPixY))
        wmin = min(wmin, poly(0))
        wmax = max(wmax, poly(NumPixY-1))
        print(i, wmin, wmax, poly(0), poly(NumPixY-1))

    #- Round down/up to nearest Angstrom
    wmin = int(wmin)
    wmax = int(wmax+1)

    #- Min and max of spot/fiber positions on the slit head
    pmin = min(spotpos[0], fiberpos[0])
    pmax = max(spotpos[-1], fiberpos[-1])

    #-------------------------------------------------------------------------
    #- For slices in wavelength, fit y vs. slit position and sample at
    #- fiberpos spoty[np, nw]

    ydeg = 7
    y_vs_w = N.zeros( (nspec, nw) )
    for i in range(nw):
        poly = Legendre.fit(spotpos, spoty[:,i], deg=ydeg, domain=(pmin, pmax))
        y_vs_w[:,i] = poly(fiberpos)

    #- For each fiber, fit y vs. wavelength and save coefficients
    #- Also calculate min/max wavelengths seen by every fiber

    wmin_all = 0
    wmax_all = 1e8
    ww = N.arange(wmin, wmax)

    ycoeff = N.zeros( (nspec, ydeg+1) )
    for i in range(nspec):
        poly = Legendre.fit(wavelength, y_vs_w[i], deg=ydeg, domain=(wmin,wmax))
        ycoeff[i] = poly.coef

        wmin_all = max(wmin_all, N.interp(0, poly(ww), ww))
        wmax_all = min(wmax_all, N.interp(NumPixY-1, poly(ww), ww))

    #- Round up/down to integer wavelengths
    wmin_all = int(wmin_all)
    wmax_all = int(wmax_all+1)

    #-------------------------------------------------------------------------
    #- for a slice in wavelength, fit x vs. slit position
    x_vs_p = N.zeros( (nw, len(fiberpos)) )
    for i in range(nw):
        poly = Legendre.fit(spotpos, spotx[:,i], deg=7, domain=(pmin, pmax))
        x_vs_p[i] = poly(fiberpos)
        assert N.max( N.abs(spotx[:,i] - poly(spotpos)) ) < 0.01

    xdeg = 7
    xcoeff = N.zeros( (nspec, xdeg+1) )
    for i in range(nspec):
        poly = Legendre.fit(wavelength, x_vs_p[:, i], deg=xdeg, domain=(wmin, wmax))
        xcoeff[i,:] = poly.coef
        assert N.max( N.abs(x_vs_p[:,i] - poly(wavelength)) ) < 0.01

    #-------------------------------------------------------------------------
    #- Write to fits file
    print("Writing", opts.outpsf)

    #- Use first spot file for representative header to pass keywords through
    hdr = fitsio.read_header(spotfiles[0])
    hdr.delete('WAVE')
    hdr.delete('FIBER')
    hdr.add_record({"name":"PSFTYPE",   "value":"SPOTGRID",     "comment":"Grid of simulated PSF spots"})
    hdr.add_record({"name":"NPIX_X",    "value":NumPixX,        "comment":"Number of CCD pixels in X direction"})
    hdr.add_record({"name":"NPIX_Y",    "value":NumPixY,        "comment":"Number of CCD pixels in Y direction"})
    hdr.add_record({"name":"NSPEC",     "value":nspec,          "comment":"Number of spectra"})
    hdr.add_record({"name":"NWAVE",     "value":nw,             "comment":"Number of wavelength samples"})
    hdr.add_record({"name":"CCDPIXSZ",  "value":CcdPixelSize,   "comment":"CCD pixel size [mm]"})
    hdr.add_record({"name":"DFIBER",    "value":FiberSpacing,   "comment":"Center-to-center pitch of fibers on slit [mm]"})
    hdr.add_record({"name":"DGROUP",    "value":GroupSpacing,   "comment":"Spacing between fiber groups on slit [mm]"})
    hdr.add_record({"name":"NGROUPS",   "value":GroupsPerCcd,   "comment":"Number of fiber groups per slit"})
    hdr.add_record({"name":"NFIBGRP",   "value":FibersPerGroup, "comment":"Number of fibers per group"})
    hdr.add_record({"name":"WAVEMIN",   "value":wmin,           "comment":"Min wavelength for Legendre domain [-1,1]"})
    hdr.add_record({"name":"WAVEMAX",   "value":wmax,           "comment":"Max wavelength for Legendre domain [-1,1]"})
    hdr.add_record({"name":"WMIN_ALL",  "value":wmin_all,       "comment":"Min wavelength seen by all spectra [Ang]"})
    hdr.add_record({"name":"WMAX_ALL",  "value":wmax_all,       "comment":"Max wavelength seen by all spectra [Ang]"})

    fitsio.write(opts.outpsf, xcoeff, extname='XCOEFF', header=hdr, clobber=True)

    wavehdr = list()
    wavehdr.append(dict(name='WAVEMIN', value=wmin, comment='Min wavelength on the CCD [Ang]'))
    wavehdr.append(dict(name='WAVEMAX', value=wmax, comment='Max wavelength on the CCD [Ang]'))
    wavehdr.append(dict(name='WMIN_ALL', value=wmin_all, comment='Min wavelength seen by all spectra [Ang]'))
    wavehdr.append(dict(name='WMAX_ALL', value=wmax_all, comment='Max wavelength seen by all spectra [Ang]'))
    fitsio.write(opts.outpsf, ycoeff, extname='YCOEFF', header=wavehdr)

    # fitsio.write(opts.outpsf, Y, extname='Y')
    # fitsio.write(opts.outpsf, W, extname='WAVELENGTH')

    fitsio.write(opts.outpsf, spots, extname='SPOTS')
    fitsio.write(opts.outpsf, spotx, extname='SPOTX')
    fitsio.write(opts.outpsf, spoty, extname='SPOTY')
    fitsio.write(opts.outpsf, fiberpos, extname='FIBERPOS')
    fitsio.write(opts.outpsf, spotpos, extname='SPOTPOS')
    fitsio.write(opts.outpsf, wavelength, extname='SPOTWAVE')

    #- Add pre-computed throughput to PSF if requested
    #- Removing; this could just lead to inconsistencies
    # if opts.throughput:
    #     header = fitsio.read_header(opts.throughput, 'THROUGHPUT')
    #     data = fitsio.read(opts.throughput, 'THROUGHPUT')
    #     fitsio.write(opts.outpsf, data, header=header, extname='THROUGHPUT')

    #--- DEBUG ---
    if opts.debug:
        import pylab as P
        P.ion()
        import IPython
        IPython.embed()
    #--- DEBUG ---
    return 0
#
#
#
exit(main())
