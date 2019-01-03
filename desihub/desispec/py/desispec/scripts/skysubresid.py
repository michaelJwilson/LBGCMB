# Script for generating plots on SkySub residuals
from __future__ import absolute_import, division

from desiutil.log import get_logger
import argparse
import numpy as np

from  desispec import _version as desis_v

def parse(options=None):
    parser = argparse.ArgumentParser(description="Generate QA on Sky Subtraction residuals [v{:s}]".format(desis_v.__offline_qa_version__))
    parser.add_argument('--reduxdir', type = str, default = None, metavar = 'PATH',
                        help = 'Override default path ($DESI_SPECTRO_REDUX/$SPECPROD) to processed data.')
    parser.add_argument('--expid', type=int, help='Generate exposure plot on given exposure')
    parser.add_argument('--night', type=str, help='Generate night plot on given night')
    parser.add_argument('--channels', type=str, help='List of channels to include')
    parser.add_argument('--prod', default=False, action="store_true", help="Results for full production run")
    parser.add_argument('--gauss', default=False, action="store_true", help="Expore Gaussianity for full production run")
    parser.add_argument('--nights', type=str, help='List of nights to limit prod plots')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args



def main(args) :
    # imports
    import glob
    from desispec.io import findfile
    from desispec.io import get_exposures
    from desispec.io import get_files, get_nights
    from desispec.io import read_frame
    from desispec.io import get_reduced_frames
    from desispec.io.sky import read_sky
    from desispec.io import specprod_root
    from desispec.qa import utils as qa_utils
    import copy
    import pdb

    # Log
    log=get_logger()
    log.info("starting")

    # Path
    if args.reduxdir is not None:
        specprod_dir = args.reduxdir
    else:
        specprod_dir = specprod_root()


    # Channels
    if args.channels is not None:
        channels = [iarg for iarg in args.channels.split(',')]
    else:
        channels = ['b','r','z']

    # Sky dict
    sky_dict = dict(wave=[], skyflux=[], res=[], count=0)
    channel_dict = dict(b=copy.deepcopy(sky_dict),
                        r=copy.deepcopy(sky_dict),
                        z=copy.deepcopy(sky_dict),
                        )

    # Exposure plot?
    if args.expid is not None:
        # Nights
        path_nights = glob.glob(specprod_dir+'/exposures/*')
        if nights is None:
            nights = get_nights()
        nights.sort()
        # Find the exposure
        for night in nights:
            if args.expid in get_exposures(night, specprod_dir=specprod_dir):
                frames_dict = get_files(filetype=str('cframe'), night=night,
                                    expid=args.expid, specprod_dir=specprod_dir)
                # Loop on channel
                #for channel in ['b','r','z']:
                for channel in ['z']:
                    channel_dict[channel]['cameras'] = []
                    for camera, cframe_fil in frames_dict.items():
                        if channel in camera:
                            sky_file = findfile(str('sky'), night=night, camera=camera,
                                expid=args.expid, specprod_dir=specprod_dir)
                            wave, flux, res, _ = qa_utils.get_skyres(cframe_fil)
                            # Append
                            channel_dict[channel]['wave'].append(wave)
                            channel_dict[channel]['skyflux'].append(np.log10(np.maximum(flux,1e-1)))
                            channel_dict[channel]['res'].append(res)
                            channel_dict[channel]['cameras'].append(camera)
                            channel_dict[channel]['count'] += 1
                    if channel_dict[channel]['count'] > 0:
                        from desispec.qa.qa_plots import skysub_resid_series  # Hidden to help with debugging
                        skysub_resid_series(channel_dict[channel], 'wave',
                             outfile='QA_skyresid_wave_expid_{:d}{:s}.png'.format(args.expid, channel))
                        skysub_resid_series(channel_dict[channel], 'flux',
                                            outfile='QA_skyresid_flux_expid_{:d}{:s}.png'.format(args.expid, channel))
        return


    # Nights
    if args.nights is not None:
        nights = [iarg for iarg in args.nights.split(',')]
    else:
        nights = None

    # Full Prod Plot?
    if args.prod:
        from desispec.qa.qa_plots import skysub_resid_dual
        # Loop on channel
        for channel in channels:
            cframes = get_reduced_frames(nights=nights, channels=[channel])
            if len(cframes) > 0:
                log.info("Loading sky residuals for {:d} cframes".format(len(cframes)))
                sky_wave, sky_flux, sky_res, _ = qa_utils.get_skyres(cframes)
                # Plot
                outfile='QA/skyresid_prod_dual_{:s}.png'.format(channel)
                log.info("Plotting to {:s}".format(outfile))
                skysub_resid_dual(sky_wave, sky_flux, sky_res, outfile=outfile)
        return

    # Full Prod Plot?
    if args.gauss:
        from desispec.qa.qa_plots import skysub_gauss
        # Loop on channel
        for channel in channels:
            cframes = get_reduced_frames(nights=nights, channels=[channel])
            if len(cframes) > 0:
                # Cut down for debugging
                #cframes = [cframes[ii] for ii in range(15)]
                #
                log.info("Loading sky residuals for {:d} cframes".format(len(cframes)))
                sky_wave, sky_flux, sky_res, sky_ivar = qa_utils.get_skyres(cframes)
                # Plot
                log.info("Plotting..")
                skysub_gauss(sky_wave, sky_flux, sky_res, sky_ivar,
                                  outfile='skyresid_prod_gauss_{:s}.png'.format(channel))
        return
