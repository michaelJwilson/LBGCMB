"""
Run integration tests from pixsim through redshifts

python -m desispec.test.integration_test
"""
from __future__ import absolute_import, print_function
import sys
import os
import random
import time
import subprocess as sp
import glob

import numpy as np
from astropy.io import fits

from desispec.util import runcmd
import desispec.pipeline as pipe
import desispec.io as io
import desiutil.log as logging

#- prevent nose from trying to run this test since it takes too long
__test__ = False


def check_env():
    """
    Check required environment variables; raise RuntimeException if missing
    """
    log = logging.get_logger()
    #- template locations
    missing_env = False
    if 'DESI_BASIS_TEMPLATES' not in os.environ:
        log.warning('missing $DESI_BASIS_TEMPLATES needed for simulating spectra'.format(name))
        missing_env = True

    if not os.path.isdir(os.getenv('DESI_BASIS_TEMPLATES')):
        log.warning('missing $DESI_BASIS_TEMPLATES directory')
        log.warning('e.g. see NERSC:/project/projectdirs/desi/spectro/templates/basis_templates/v2.2')
        missing_env = True

    for name in (
        'DESI_SPECTRO_SIM', 'DESI_SPECTRO_REDUX', 'PIXPROD', 'SPECPROD', 'DESIMODEL'):
        if name not in os.environ:
            log.warning("missing ${0}".format(name))
            missing_env = True

    if missing_env:
        log.warning("Why are these needed?")
        log.warning("    Simulations written to $DESI_SPECTRO_SIM/$PIXPROD/")
        log.warning("    Raw data read from $DESI_SPECTRO_DATA/")
        log.warning("    Spectro pipeline output written to $DESI_SPECTRO_REDUX/$SPECPROD/")
        log.warning("    Templates are read from $DESI_BASIS_TEMPLATES")

    #- Wait until end to raise exception so that we report everything that
    #- is missing before actually failing
    if missing_env:
        log.critical("missing env vars; exiting without running pipeline")
        sys.exit(1)

    #- Override $DESI_SPECTRO_DATA to match $DESI_SPECTRO_SIM/$PIXPROD
    os.environ['DESI_SPECTRO_DATA'] = os.path.join(os.getenv('DESI_SPECTRO_SIM'), os.getenv('PIXPROD'))


# Simulate raw data

def sim(night, nspec=5, clobber=False):
    """
    Simulate data as part of the integration test.

    Args:
        night (str): YEARMMDD
        nspec (int, optional): number of spectra to include
        clobber (bool, optional): rerun steps even if outputs already exist

    Raises:
        RuntimeError if any script fails
    """
    log = logging.get_logger()

    # Create input fibermaps, spectra, and pixel-level raw data

    for expid, program in zip([0,1,2], ['flat', 'arc', 'dark']):
        cmd = "newexp-random --program {program} --nspec {nspec} --night {night} --expid {expid}".format(
            expid=expid, program=program, nspec=nspec, night=night)
        fibermap = io.findfile('fibermap', night, expid)
        simspec = '{}/simspec-{:08d}.fits'.format(os.path.dirname(fibermap), expid)
        inputs = []
        outputs = [fibermap, simspec]
        if runcmd(cmd, inputs=inputs, outputs=outputs, clobber=clobber) != 0:
            raise RuntimeError('newexp-random failed for {} exposure {}'.format(program, expid))

        cmd = "pixsim --preproc --nspec {nspec} --night {night} --expid {expid}".format(expid=expid, nspec=nspec, night=night)
        inputs = [fibermap, simspec]
        outputs = list()
        outputs.append(fibermap.replace('fibermap-', 'simpix-'))
        for camera in ['b0', 'r0', 'z0']:
            pixfile = io.findfile('pix', night, expid, camera)
            outputs.append(pixfile)
            #outputs.append(os.path.join(os.path.dirname(pixfile), os.path.basename(pixfile).replace('pix-', 'simpix-')))
        if runcmd(cmd, inputs=inputs, outputs=outputs, clobber=clobber) != 0:
            raise RuntimeError('pixsim failed for {} exposure {}'.format(program, expid))

    return


def integration_test(night=None, nspec=5, clobber=False):
    """Run an integration test from raw data simulations through redshifts

    Args:
        night (str, optional): YEARMMDD, defaults to current night
        nspec (int, optional): number of spectra to include
        clobber (bool, optional): rerun steps even if outputs already exist

    Raises:
        RuntimeError if any script fails

    """
    log = logging.get_logger()
    log.setLevel(logging.DEBUG)

    # YEARMMDD string, rolls over at noon not midnight
    # TODO: fix usage of night to be something other than today
    if night is None:
        #night = time.strftime('%Y%m%d', time.localtime(time.time()-12*3600))
        night = "20160726"

    # check for required environment variables
    check_env()

    # simulate inputs
    sim(night, nspec=nspec, clobber=clobber)

    # create production

    # FIXME:  someday run PSF estimation too...
    ### com = "desi_pipe --spectrographs 0 --fakeboot --fakepsf"
    com = "desi_pipe --spectrographs 0 --fakeboot --fakepsf"
    sp.check_call(com, shell=True)

    # raw and production locations

    rawdir = os.path.abspath(io.rawdata_root())
    proddir = os.path.abspath(io.specprod_root())

    # Modify options file to restrict the spectral range

    optpath = os.path.join(proddir, "run", "options.yaml")
    opts = pipe.yaml_read(optpath)
    opts['extract']['specmin'] = 0
    opts['extract']['nspec'] = nspec
    pipe.yaml_write(optpath, opts)

    # run the generated shell scripts

    # FIXME:  someday run PSF estimation too...

    # print("Running bootcalib script...")
    # com = os.path.join(proddir, "run", "scripts", "bootcalib_all.sh")
    # sp.check_call(["bash", com])

    # print("Running specex script...")
    # com = os.path.join(proddir, "run", "scripts", "specex_all.sh")
    # sp.check_call(["bash", com])

    # print("Running psfcombine script...")
    # com = os.path.join(proddir, "run", "scripts", "psfcombine_all.sh")
    # sp.check_call(["bash", com])

    com = os.path.join(proddir, "run", "scripts", "run_shell.sh")
    print("Running extraction through calibration: "+com)
    sp.check_call(["bash", com])

    com = os.path.join(proddir, "run", "scripts", "spectra.sh")
    print("Running spectral regrouping: "+com)
    sp.check_call(["bash", com])

    com = os.path.join(proddir, "run", "scripts", "redshift.sh")
    print("Running redshift script "+com)
    sp.check_call(["bash", com])

    # #-----
    # #- Did it work?
    # #- (this combination of fibermap, simspec, and zbest is a pain)
    expid = 2
    fmfile = io.findfile('fibermap', night=night, expid=expid)
    fibermap = io.read_fibermap(fmfile)
    simdir = os.path.dirname(fmfile)
    simspec = '{}/simspec-{:08d}.fits'.format(simdir, expid)
    siminfo = fits.getdata(simspec, 'TRUTH')

    from desimodel.footprint import radec2pix
    nside=64
    pixels = np.unique(radec2pix(nside, fibermap['RA_TARGET'], fibermap['DEC_TARGET']))

    print()
    print("--------------------------------------------------")
    print("Pixel     True  z        ->  Class  z        zwarn")
    # print("3338p190  SKY   0.00000  ->  QSO    1.60853   12   - ok")
    for pix in pixels:
        zbest = io.read_zbest(io.findfile('zbest', groupname=pix))
        for i in range(len(zbest.z)):
            objtype = zbest.spectype[i]
            z, zwarn = zbest.z[i], zbest.zwarn[i]

            j = np.where(fibermap['TARGETID'] == zbest.targetid[i])[0][0]
            truetype = siminfo['OBJTYPE'][j]
            oiiflux = siminfo['OIIFLUX'][j]
            truez = siminfo['REDSHIFT'][j]
            dv = 3e5*(z-truez)/(1+truez)
            if truetype == 'SKY' and zwarn > 0:
                status = 'ok'
            elif truetype == 'ELG' and zwarn > 0 and oiiflux < 8e-17:
                status = 'ok ([OII] flux {:.2g})'.format(oiiflux)
            elif zwarn == 0:
                if truetype == 'LRG' and objtype == 'GALAXY' and abs(dv) < 150:
                    status = 'ok'
                elif truetype == 'ELG' and objtype == 'GALAXY':
                    if abs(dv) < 150:
                        status = ok
                    elif oiiflux < 8e-17:
                        status = 'ok ([OII] flux {:.2g})'.format(oiiflux)
                    else:
                        status = 'OOPS ([OII] flux {:.2g})'.format(oiiflux)
                elif truetype == 'QSO' and objtype == 'QSO' and abs(dv) < 750:
                    status = 'ok'
                elif truetype in ('STD', 'FSTD') and objtype == 'STAR':
                    status = 'ok'
                else:
                    status = 'OOPS'
            else:
                status = 'OOPS'
            print('{0:<8d}  {1:4s} {2:8.5f}  -> {3:5s} {4:8.5f} {5:4d}  - {6}'.format(
                pix, truetype, truez, objtype, z, zwarn, status))

    print("--------------------------------------------------")




if __name__ == '__main__':
    integration_test()
