from    __future__           import  absolute_import, division, print_function

import  json
import  os, sys
import  warnings
import  numpy                as      np

import  rrio
import  desispec.io

from    collections          import  OrderedDict 
from    zfind                import  zfind
from    desispec.resolution  import  Resolution
from    dataobj              import  Target, MultiprocessingSharedSpectrum, SimpleSpectrum, MPISharedTargets
from    astropy.table        import  Column


if sys.version_info[0] > 2:
  basestring = str

def write_zbest(outfile, zbest, overwrite=True):
    '''
    Write zbest Table to outfile.

    Adds blank BRICKNAME and SUBTYPE columns if needed
    Adds zbest.meta['EXTNAME'] = 'ZBEST'
    '''
    ntargets = len(zbest)

    if 'BRICKNAME' not in zbest.colnames:
        zbest['BRICKNAME'] = np.zeros(ntargets, dtype='S8')

    if 'SUBTYPE' not in zbest.colnames:
        zbest['SUBTYPE']   = np.zeros(ntargets, dtype='S8')

    zbest.meta['EXTNAME']  = 'ZBEST'

    zbest.write(outfile, overwrite=overwrite)

def read_spectra(spectrafiles, targetids=None, spectrum_class=MultiprocessingSharedSpectrum):
    '''
    Read targets from a list of spectra files.
   
    Args:
        spectrafiles : list of input spectra files, or string glob to match

    Options:
        targetids :       List of target ids. If set, only those target spectra will be read.
        spectrum_class :  
        The spectrum_class argument is needed to use the same read_spectra
        routine for the two parallelism approaches for redrock. The multiprocessing version
        uses a shared memory at the spectrum class initialization, see
        the redrock.dataobj.MultiprocessingSharedSpectrum class, whereas the MPI version
        implements the shared memory after all spectra have been read by the root process,
        and so the MPI version used another more simple spectrum class (see redrock.dataobj.SimpleSpectrum).
    
    Returns tuple of (targets, meta) where
        targets is a list of Target objects and
        meta is a Table of metadata (currently only BRICKNAME)
    '''

    assert len(spectrafiles) > 0

    input_spectra   = list()
    input_targetids =  set()

    print("(REDROCK) Loading spectra files:  ", spectrafiles)

    ##  Ignore warnings about zdc2 bricks lacking bricknames in header
    for i, infile in enumerate(spectrafiles):
        sp = desispec.io.read_spectra(infile)

        if hasattr(sp, 'fmap'):
            sp.fibermap = sp.fmap 

        input_spectra.append(sp)
        
        ## print("\n\nWarning:  Overwriting TARGETIDS from %s with enumeration: %d \n\n" % (infile, i))
        ## sp.fibermap['TARGETID'] = Column(data=[i], name='TARGETID')

        input_targetids.update(sp.fibermap['TARGETID'])
        
    if targetids is None:
        targetids = input_targetids

        print("Target ids of input specfile:  ", targetids)

    targets    = list()
    bricknames = list()

    for targetid in targetids:
        spectra = list()

        for sp in input_spectra:
            ii  = (sp.fibermap['TARGETID'] == targetid)

            if np.count_nonzero(ii) == 0:
                continue

            if 'BRICKNAME' in sp.fibermap.dtype.names:
                brickname = sp.fibermap['BRICKNAME'][ii][0]

            else:
                brickname = 'unknown'

            for x in sp.bands:  ## typically 'b', 'r', 'z'
                wave  = sp.wave[x]                
                flux  = sp.flux[x][ii]
                ivar  = sp.ivar[x][ii]
                Rdata = sp.resolution_data[x][ii]

                for i in range(flux.shape[0]):
                    if np.all(flux[i] == 0):
                        continue

                    if np.all(ivar[i] == 0):
                        continue

                    R = Resolution(Rdata[i])

                    spectra.append(spectrum_class(wave, flux[i], ivar[i], R))

        bricknames.append(brickname)

        if len(spectra) > 0:
            targets.append(Target(targetid, spectra))

        else:
            print('ERROR: Target {} on {} has no good spectra'.format(targetid, os.path.basename(brickfiles[0])))

    ##  Create a metadata table in case we might want to add other columns in the future
    assert len(bricknames) == len(targets)

    meta                   = OrderedDict()

    meta['BRICKNAME']      = bricknames
    meta['SPECTRAFILES']   = spectrafiles

    return  targets, meta

def read_bricks(brickfiles, trueflux=False, targetids=None, spectrum_class=SimpleSpectrum):
    '''
    Read targets from a list of brickfiles

    Args:
        brickfiles : list of input brick files, or string glob to match
    
    Options:
        targetids : list of target ids. If set, only those target spectra will be read.
        spectrum_class :  The spectrum_class argument is needed to use the same read_spectra
                          routine for the two parallelism approaches for redrock. The multiprocessing version
                          uses a shared memory at the spectrum class initialization, see
                          the redrock.dataobj.MultiprocessingSharedSpectrum class, whereas the MPI version
                          implements the shared memory after all spectra have been read by the root process,
                          and so the MPI version used another more simple spectrum class (see redrock.dataobj.SimpleSpectrum).
    
    Returns list of Target objects

    Note: these don't actually have to be bricks anymore; they are read via
        desispec.io.read_frame()
    '''
    if isinstance(brickfiles, basestring):
        import glob

        brickfiles = glob.glob(brickfiles)

    assert len(brickfiles) > 0

    bricks = list()
    brick_targetids = set()

    ## Ignore warnings about zdc2 bricks lacking bricknames in header
    for infile in brickfiles:
        b = desispec.io.read_frame(infile)
        bricks.append(b)
        brick_targetids.update(b.fibermap['TARGETID'])

    if targetids is None:
        targetids = brick_targetids

    targets = list()
    bricknames = list()

    for targetid in targetids:
        spectra = list()

        for brick in bricks:
            wave = brick.wave
            ii = (brick.fibermap['TARGETID'] == targetid)

            if np.count_nonzero(ii) == 0:
                continue

            if 'BRICKNAME' in brick.fibermap.dtype.names:
                brickname = brick.fibermap['BRICKNAME'][ii][0]

            else:
                brickname = 'unknown'
            flux = brick.flux[ii]
            ivar = brick.ivar[ii]
            Rdata = brick.resolution_data[ii]

            ## work around desispec.io.Brick returning 32-bit non-native endian
            # flux = flux.astype(float)
            # ivar = ivar.astype(float)
            # Rdata = Rdata.astype(float)

            for i in range(flux.shape[0]):
                if np.all(flux[i] == 0):
                    print('WARNING: Skipping spectrum {} of target {} on brick {} with flux=0'.format(i, targetid, brick.brickname))
                    continue

                if np.all(ivar[i] == 0):
                    print('WARNING: Skipping spectrum {} of target {} on brick {} with ivar=0'.format(i, targetid, brick.brickname))
                    continue

                R = Resolution(Rdata[i])
                spectra.append(spectrum_class(wave, flux[i], ivar[i], R))

        if len(spectra) > 0:
            bricknames.append(brickname)
            targets.append(Target(targetid, spectra))

        else:
            print('ERROR: Target {} on {} has no good spectra'.format(targetid, os.path.basename(brickfiles[0])))

    ## Create a metadata table in case we might want to add other columns
    ## in the future
    assert len(bricknames) == len(targets)

    dtype = [('BRICKNAME', 'S8'),]

    meta = np.zeros(len(bricknames), dtype=dtype)

    meta['BRICKNAME'] = bricknames

    return targets, meta


if  __name__ == "__main__":
    import  time 
    import  optparse    

    from    astropy.io  import fits


    print("\n\nWelcome to DESI redrock.\n\n")
        
    start_time       = time.time()    
    
    ncpu             = 1
    survey           = 'desi'  ## 'pfs'

    ## redshifts     = 3.5 + np.linspace(0., 2.50, 25)
    ## exposures     = 900 * np.arange(1, 22, 1)

    redshifts        = np.array([2.0])
    exposures        = 60. * 15. * np.arange(1, 2, 1)

    ##  type         = 'BC03'
    ##  subpath      = ''

    ##  Shapley quantile number. 
    type             = 'Shapley'
    quantile         = 3
    subpath          = type + '/Q%d/' % quantile

    t0               = time.time()

    templates        = rrio.read_templates()         
    
    for redshift in redshifts:
      for exposure in exposures:
        ## infiles       = [os.environ['LBGCMB'] + '/quickspectra/dat/BC03/%s/spec-BC03-z%.1lf_exp%d.fits' % (survey, redshift, exposure)]
        infiles          = [os.environ['LBGCMB'] + '/quickspectra/dat/%s/spec-Shapley-Q3_199-z%.1lf_exp%d.fits' % (subpath, redshift, exposure)] 

        output           = infiles[0].split('/')[-1].split('.fits')[0]
        
        try:
          targets, meta  =  read_spectra(infiles) ## spectrum_class=SimpleSpectrum
        
        except RuntimeError:
          targets, meta  =  read_bricks(infiles) ## spectrum_class=SimpleSpectrum
        
        dt               = time.time() - t0
    
        zscan, zfit      = zfind(targets, templates, ncpu=ncpu)
        
        rrio.write_zscan('./dat/rrh5/%s/%s/%s.h5' % (survey, subpath, output), zscan, zfit, clobber=True)
    
        zbest = zfit[zfit['znum'] == 0]

        ##  Remove extra columns not needed for zbest;  'zzchi2'
        zbest.remove_columns(['zz']) 

        ##  Change to upper case like DESI
        for colname in zbest.colnames:
          if colname.islower():
            zbest.rename_column(colname, colname.upper())

        ##  Add brickname column
        zbest['BRICKNAME'] = meta['BRICKNAME']

        write_zbest('./dat/zbest/%s/%s/%s.fits' % (survey, subpath, output), zbest)
    
        run_time = time.time() - start_time
    
        json.dump(meta, open('./dat/meta/%s/%s/%s_meta.json' % (survey, subpath, output), 'w'))    
            
    print("\n\nDone.\n\n")
