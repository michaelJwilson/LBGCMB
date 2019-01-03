#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compute the Lyman-alpha S/N for the mean QSO at different g-mag, z,
as a function of wavelength.
It is the python version of the IDL code desi_quicklya.pro by D. Schlegel (LBL),
and it is heavily based on quicksim.py by D. Kirkby (UC Irvine)

Use desi_quicklya.py --help for instructions on running this program.

The DESIMODEL environment variable should be set to the root of your
desimodel package and PYTHONPATH should include $DESIMODEL/py/, e.g.

export DESIMODEL=`pwd`
export PYTHONPATH=$DESIMODEL/py/:$PYTHONPATH

Created 06-May-2015 by Andreu Font-Ribera (afont@lbl.gov)
"""

import argparse
import os
import os.path
import numpy as np

import desimodel.simulate as sim

def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v','--verbose', action = 'store_true',
        help = 'provide verbose output on progress')
    parser.add_argument('--ab-magnitude', type = str, default = "g=23.0",
                    help = 'max magnitude to compute, e.g. g=22.0 or r=21.5')
    parser.add_argument('--exptime', type = float, default = 4000,
        help = 'overrides exposure time specified in the parameter file (secs)')
    args = parser.parse_args()

    # We require that the DESIMODEL environment variable is set.
    if 'DESIMODEL' not in os.environ:
        raise RuntimeError('The environment variable DESIMODEL must be set.')

    # Load the source spectrum to use.
    infile=os.environ['DESIMODEL']+'/data/spectra/spec-lya.dat'
    if not os.path.isfile(infile):
        print('Could not find file: %s' % infile)
        return -1

    # these were options in quicksim.py, but we keep them fix here (as in IDL)
    sky='dark'
    airmass=1.0
    min_wavelength=3500.25
    max_wavelength=9999.7
    wavelength_step=1.0
    downsampling=1.0
    # specify target type (in the IDL version this is star, they are the same)
    model='qso'

    # figure out magnitude
    try:
        band = args.ab_magnitude[0]
        abmag = float(args.ab_magnitude[2:])
        assert band in 'ugriz' and args.ab_magnitude[1] == '='
    except(AssertionError,ValueError):
        print('Invalid ab-magnitude parameter. '
                                +'Valid syntax is, e.g. g=22.0 or r=21.5.')
        return -1

    # will generate a file for each g magnitude
    min_m=19.25
    dm=0.5
    Nm=np.ceil((abmag-min_m)/dm)
    mags = np.linspace(min_m, min_m+Nm*dm, Nm+1)
    if args.verbose: print('mags', mags)

    # will loop over these redshifts, compute SN for each and collect them
    # this is the same grid in the infile file, can't be changed here!
    zqs = np.linspace(2.0, 4.75, 12)

    for mag in mags:
        if args.verbose: print('generate file for %s = %f' % (band,mag) )

        # for each wavelength, we will store a SN per redshift
        collect_results = []

        for i,zq in enumerate(zqs):
            if args.verbose: print('compute SN for zq = %f' % zq)
            # wavelength colums is always the first one
            wave_col=0
            # after wavelength we have the flux template for the different z
            flux_col=i+1
            srcSpectrum = sim.SpectralFluxDensity.loadFromTextFile(infile,
                            wavelengthColumn=wave_col,valuesColumn=flux_col,
                            extrapolatedValue=(0.))
            # Rescale the spectrum to the requested magnitude
            if args.verbose:
                print('Rescaling %s-band magnitude to %f' % (band,mag) )
            srcSpectrum = srcSpectrum.createRescaled(band,mag)

            # Calculate and print different AB magnitudes of the source spectrum
            # (not used anywhere else in the code)
            if args.verbose:
                ms = srcSpectrum.getABMagnitudes()
                specSummary = 'zq = %.2f quasar, ' % zq
                for b in 'ugriz':
                    # Check for a valid AB magnitude in this band.
                    if ms[b] is not None: specSummary += ' %s=%.2f' % (b,ms[b])
                print(specSummary)

            # Create the default atmosphere for the requested sky conditions.
            atmosphere = sim.Atmosphere(skyConditions=sky,
                                            basePath=os.environ['DESIMODEL'])

            # Create a quick simulator using the default instrument model.
            qsim = sim.Quick(atmosphere=atmosphere,
                                            basePath=os.environ['DESIMODEL'])

            # Initialize the simulation wavelength grid to use.
            qsim.setWavelengthGrid(min_wavelength,max_wavelength,
                                                                wavelength_step)

            # Perform a quick simulation of the observed spectrum.
            if args.verbose: print('Running quick simulation.')
            results = qsim.simulate(sourceType=model,
                        sourceSpectrum=srcSpectrum,airmass=airmass,
                        expTime=args.exptime,downsampling=downsampling)

            # store values for each wavelength
            for j,row in enumerate(results):
                if i==0:
                    new_row = {'wave': row.wave, 'snr_'+str(zq): row.snrtot}
                    collect_results.append(new_row)
                else:
                    collect_results[j]['snr_'+str(zq)] = row.snrtot

            # Calculate the median total SNR in bins with some observed flux.
            medianSNR = np.median(results[results.obsflux > 0].snrtot)
            # Calculate the total SNR^2 for the combined cameras.
            totalSNR2 = np.sum(results.snrtot**2)
            # Print a summary of SNR statistics.
            if args.verbose:
                snrSummary = 'Median S/N = %.3f, ' % medianSNR
                snrSummary += ' Total (S/N)^2 = %.1f' % totalSNR2
                print(snrSummary)

        # Save the results to file
        fname = os.environ['DESIMODEL']+'/data/spectra/sn-spec-lya-'
        fname += band+str(mag)+'-t'+str(int(args.exptime))+'.dat'
        if args.verbose: print('Saving results to %s' % fname)
        # Try opening the requested output file.
        with open(fname,'w') as out:
            print('# Lyman-alpha forest S/N per Ang for mean quasar with '
                                                      +'mean forest',file=out)
            print('# INFILE=',infile,file=out)
            print('# BAND=',band,file=out)
            print('# MAG=',mag,file=out)
            print('# EXPTIME=',qsim.expTime,file=out)
            print('#',file=out)
            z_header = '# Wave'
            for zq in zqs: z_header += ' SN(z='+str(zq)+')'
            print(z_header,file=out)

            for row in collect_results:
                toprint = str('%9.2f' % row['wave'])
                for zq in zqs:
                    toprint += ' %9.2f' % row['snr_'+str(zq)]
                print(toprint, file=out)

if __name__ == '__main__':
    main()
