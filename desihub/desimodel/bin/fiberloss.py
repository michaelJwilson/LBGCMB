#!/usr/bin/env python
#
# python version of fiberloss.pro
# J. Guy 2014/06/23
#
# adapted to include telescope blur and updated wavelength dependence
# S. Bailey 2014/06/26
#
# Compute the fiber input losses for the benchmark spectra cases
# Fiber diameter 1.5 arcsec
# Seeing is 1.1 arcsec FWHM with Moffat beta=3.5, and scale with
# wavelength as (6355./wave)^0.2
# Lateral offset is 0.2 arcsec from fiber center
# Cases are:
#   * Point source
#   * ELG with half-light radius 0.45 arcsec
#   * LRG with half-light radius 1.0 arcsec
#
# TO DO
# - Refactor into module function + wrapper script
# - Extract hard coded parameters into a fiberloss.yaml file

import sys, os
from math import *
import numpy
from scipy.signal import fftconvolve
from scipy.interpolate import InterpolatedUnivariateSpline
import yaml

#- Check options before doing anything
import argparse
parser = argparse.ArgumentParser(prog=sys.argv[0])
parser.add_argument("-o", "--outdir", action='store', metavar='DIR', help="write fiberloss-*.dat files to DIR.")
parser.add_argument("-t", "--test", action="store_true", help="test convolution code")
opts = parser.parse_args()


#- Read desi.yaml to get typical fiber size and plate scale
desiparams = yaml.load(open(os.getenv('DESIMODEL')+'/data/desi.yaml'))
platescale = desiparams['fibers']['diameter_um'] / desiparams['fibers']['diameter_arcsec']
fiber_radius = desiparams['fibers']['diameter_arcsec'] / 2.0

#- read DESI blur from 0347 system throughput file
blurfile = os.getenv('DESIMODEL')+'/data/inputs/throughput/DESI-0347-throughput.txt'
tmp = numpy.loadtxt(blurfile).T
blurwave = tmp[0] * 10   #- nm -> Angstroms
blur_arcsec = tmp[4] / platescale
blur = InterpolatedUnivariateSpline(blurwave, blur_arcsec)

#- Hard coded parameters; should move into a parameter file
seeing       = 1.1  # FWHM in arcsec
wave_ref     = 6355 # reference wavelength for seeing
seeing_scale = 0.2  # seeing(wave) = seeing * (wave_ref/wave)^seeing_scale
beta         = 3.5  # Moffat profile parameter
offset       = 0.20 # lateral error in arcsec
r_elg        = 0.45 # half-light radius for ELGs
r_lrg        = 1.0  # half-light radius for LRGs
### mayall_blur  = 0.22 # RMS blur (arcsec) from telescope included in Dey & Valdez seeing
mayall_blur = desiparams['jacoby_seeing']

#- size and resolution of image to simulate
npix         = 501   # number of pixels across image
pixscale     = 0.025 # arcsec/pix

#- Remove mayall_blur from seeing
seeing_sigma = numpy.sqrt( (seeing/2.35482)**2 - mayall_blur**2 )
site_seeing = seeing_sigma * 2.35482

#- Coordinates and mask for offset fiber on the image
coord  = (numpy.arange((npix))-npix/2)*pixscale # coord along one axis, in arcsec
x      = coord.reshape(npix,1)
y      = coord.reshape(1,npix)
x      = x+0*y # now a 2D array of the x coordinates
y      = y+0*x # now a 2D array of the y coordinates

radius                   = numpy.sqrt(x**2+y**2) # distance from center of object
radius_from_fiber_center = numpy.sqrt((x-offset)**2+y**2) # distance from center of fiber

fiber_mask = (radius_from_fiber_center<fiber_radius) # 1=in fiber


# test of numpy.fft.rfft2 or fftconvolve, who knows what happens with another version of python
if opts.test :
    sigma=1.1
    img_seeing   = 1./(2*pi*sigma**2)*numpy.exp(-(radius**2/2/sigma**2))
    # img_seeing2  = numpy.fft.irfft2(numpy.fft.rfft2(img_seeing) * numpy.fft.rfft2(img_seeing), img_seeing.shape)
    # img_seeing2  = numpy.roll(numpy.roll(img_seeing2,npix/2+1,axis=0),npix/2+1,axis=1)*pixscale**2 # re-center and normalize
    img_seeing2  = fftconvolve(img_seeing, img_seeing, mode='same')*pixscale**2
    img_seeing3  = 1./(2*pi*2*sigma**2)*numpy.exp(-(radius**2/2/(2*sigma**2)))
    if abs(numpy.sum(img_seeing2)/numpy.sum(img_seeing3)-1)>0.01 or abs(numpy.sum(img_seeing2*fiber_mask)/numpy.sum(img_seeing3*fiber_mask)-1)>0.01:
        print("error in the convolution using numpy.fft.rfft2")
        print(numpy.sum(img_seeing)*pixscale**2,numpy.sum(img_seeing2)*pixscale**2,numpy.sum(img_seeing3)*pixscale**2)
        sys.exit(12)
    else:
        print("Convolution test passed; your scipy is OK")
        sys.exit(0)

waves      = numpy.arange(3500.,10001.,250.)
input_sky  = numpy.ones(len(waves))
input_star = numpy.zeros(len(waves))
input_elg  = numpy.zeros(len(waves))
input_lrg  = numpy.zeros(len(waves))

print("# Fiber input geometric acceptance")
print("# {} arcsec seeing at {} Angstroms, Moffat beta={}".format(seeing, wave_ref, beta))
print("# seeing wavelength dependence scales as ({}/lambda)^{}".format(wave_ref, seeing_scale))
print("# removes {} arcsec existing Mayall blur then adds DESI blur".format(mayall_blur))
print("# lateral positioner offset {} arcsec".format(offset))
print("# point = point source, e.g. a star")
print("# elg   = Exponential half-light radius {}".format(r_elg))
print("# lrg   = de Vaucouleurs half-light radius {}".format(r_lrg))
print("# ")
print("# wave   star     elg      lrg")

#- Reference ELG and LRG images
img_elg    = numpy.exp(-1.678*radius/r_elg)
img_lrg = numpy.exp(-7.67*((radius/r_lrg)**0.25 - 1))

for i, wave in enumerate(waves):
    #- Input seeing image
    fwhm       = site_seeing * (wave_ref/wave)**seeing_scale
    alpha      = 0.5 * fwhm / sqrt(2.**(1./beta) -1)
    img_seeing = (beta - 1) * (pi * alpha**2)**(-1) * (1 + radius**2/alpha**2)**(-beta)

    #- convolve with telescope blur
    img_blur   = numpy.exp(-radius**2 / (2*blur(wave)**2))
    img_blur  /= numpy.sum(img_blur)
    img_seeing = fftconvolve(img_seeing, img_blur, mode='same')

    #- Point source (star) input acceptance
    input_star[i] = numpy.sum(img_seeing*fiber_mask)/numpy.sum(img_seeing)

    # for comparison with fiberloss.pro, trim img_seeing
    # img_seeing *= ((abs(x)<0.2*npix*pixscale)&(abs(y)<0.2*npix*pixscale))

    #- Convolve ELG image with seeing+blur
    #img_elg2   = numpy.fft.irfft2(numpy.fft.rfft2(img_elg) * numpy.fft.rfft2(img_seeing), img_elg.shape)
    #img_elg2   = numpy.roll(numpy.roll(img_elg2,npix/2+1,axis=0),npix/2+1,axis=1)*pixscale**2 # re-center and normalize
    img_elg2 = fftconvolve(img_elg, img_seeing, mode='same') * pixscale**2

    #- ELG input acceptance
    input_elg[i] = numpy.sum(img_elg2*fiber_mask)/numpy.sum(img_elg2)

    #- Convolve LRG image with seeing+blur
    #img_lrg2   = numpy.fft.irfft2(numpy.fft.rfft2(img_lrg) * numpy.fft.rfft2(img_seeing), img_lrg.shape)
    #img_lrg2   = numpy.roll(numpy.roll(img_lrg2,npix/2+1,axis=0),npix/2+1,axis=1)*pixscale**2 # re-center and normalize
    img_lrg2 = fftconvolve(img_lrg, img_seeing, mode='same') * pixscale**2

    #- LRG input acceptance
    input_lrg[i] = numpy.sum(img_lrg2*fiber_mask)/numpy.sum(img_lrg2)

    print("{:7.1f} {:8.5f} {:8.5f} {:8.5f}".format(wave, input_star[i], input_elg[i], input_lrg[i]))

#- Write output files if --outdir is given
if opts.outdir is not None:
    header = """\
    # Fiber input geometric acceptance
    # {seeing} arcsec seeing with Moffat beta={beta} scales as ({wave_ref}/lambda)^{seeing_scale}
    # Calculation removes {mayall_blur} arcsec existing Mayall blur then adds DESI blur
    # Lateral positioner offset={offset} arcsec""".format(
        seeing=seeing, beta=beta, wave_ref=wave_ref, seeing_scale=seeing_scale,
        mayall_blur=mayall_blur, offset=offset)

    #- Point source / star
    fx = open(opts.outdir+'/fiberloss-star.dat', 'w')
    print(header, file=fx)
    print("#", file=fx)
    print("# Star / point-source", file=fx)
    print("#", file=fx)
    print("# Wavelength FiberAcceptance", file=fx)
    for w, t in zip(waves, input_star):
        print("{:7.1f} {:8.5f}".format(w, t), file=fx)
    fx.close()

    #- Also write out QSOs as if they are point sources
    fx = open(opts.outdir+'/fiberloss-qso.dat', 'w')
    print(header, file=fx)
    print("#", file=fx)
    print("# QSO, treating as a point-source", file=fx)
    print("#", file=fx)
    print("# Wavelength FiberAcceptance", file=fx)
    for w, t in zip(waves, input_star):
        print("{:7.1f} {:8.5f}".format(w, t), file=fx)
    fx.close()

    fx = open(opts.outdir+'/fiberloss-elg.dat', 'w')
    print(header, file=fx)
    print("#", file=fx)
    print("# ELG with exponential half-light radius={} arcsec".format(r_elg), file=fx)
    print("#", file=fx)
    print("# Wavelength FiberAcceptance", file=fx)
    for w, t in zip(waves, input_elg):
        print("{:7.1f} {:8.5f}".format(w, t), file=fx)
    fx.close()

    fx = open(opts.outdir+'/fiberloss-lrg.dat', 'w')
    print(header, file=fx)
    print("#", file=fx)
    print("# LRG with de Vaucouleurs half-light radius={} arcsec".format(r_lrg), file=fx)
    print("#", file=fx)
    print("# Wavelength FiberAcceptance", file=fx)
    for w, t in zip(waves, input_lrg):
        print("{:7.1f} {:8.5f}".format(w, t), file=fx)
    fx.close()

    fx = open(opts.outdir+'/fiberloss-sky.dat', 'w')
    print("# Sky and calibration lamp spectra do not have fiber input losses", file=fx)
    print("#", file=fx)
    print("# Wavelength FiberAcceptance", file=fx)
    for w, t in zip(waves, input_sky):
        print("{:7.1f} {:8.5f}".format(w, t), file=fx)
    fx.close()
