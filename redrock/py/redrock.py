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
from    rr                   import  write_zbest, read_spectra, read_bricks


if __name__ == '__main__':
    import  time 
    import  optparse    

    from    astropy.io       import  fits


    print("\n\nWelcome to DESI redrock.\n\n")
        
    start_time       =  time.time()    
    
    ncpu             =  1

    survey           =  'desi'
    ## survey        =  'pfs'

    ##  Sub-sample redshifts and exposures (not magnitudes).
    sample_rate      =  1

    ##  BC03 production run. 
    ##  redshifts    =  3.5 + np.linspace(0., 2.50, 25)
    ##  exposures    =  900 * np.arange(1, 22, 1)

    ##  Shapley production run. 
    redshifts        =   1.5 + np.linspace(0., 1.50, 15)[::sample_rate]
    magnitudes       =  20.0 +   np.arange(0., 6.00,  1)    
    exposures        =  60.  * 15. * np.arange(1, 30, 3)[::sample_rate]

    print('\n\nSolving for redshifts:\n'  + ''.join('%.3lf, ' % x for x in redshifts))
    print('\n\nSolving for magnitudes:\n' + ''.join('%.3lf, ' % x for x in magnitudes))
    print('\n\nSolving for exposures:\n'  + ''.join('%.3lf, ' % x for x in exposures))

    ##  type         =  'BC03' 
    ##  subpath      =  ''

    ##  Shapley quantile number. 
    type             =  'Shapley'
    Quantile         =  3
    subpath          =  type + '/all'

    t0               =  time.time()

    templates        =  rrio.read_templates()         
    nruns            =  len(exposures) * len(redshifts)

    for ii, exposure in enumerate(exposures):
      for jj, redshift in enumerate(redshifts):    
        ## infiles  = [os.environ['LBGCMB'] + '/quickspectra/dat/quickspectra/%s/%s/spec-BC03-z%.1lf_exp%d.fits' % (survey, type, redshift, exposure)]
        infiles     = [os.environ['LBGCMB'] + '/quickspectra/dat/quickspectra/%s/%s/spec-Shapley-all-z%.3lf_exp%d.fits' % (survey, subpath,\
                                                                                                                           redshift, exposure)] 
        output           =  infiles[0].split('/')[-1].split('.fits')[0]
        
        try:
          targets, meta  =  read_spectra(infiles) ## spectrum_class=SimpleSpectrum
        
        except  RuntimeError:
          targets, meta  =  read_bricks(infiles) ## spectrum_class=SimpleSpectrum
        
        dt               =  time.time() - t0
        
        zscan, zfit      =  zfind(targets, templates, ncpu=ncpu)
        

        rrio.write_zscan(os.environ['LBGCMB'] + '/redrock/dat/rrh5/%s/%s/%s.h5' % (survey, subpath, output), zscan, zfit, clobber=True)
        
        zbest = zfit[zfit['znum'] == 0]

        ##  Remove extra columns not needed for zbest;  'zzchi2'
        zbest.remove_columns(['zz']) 

        ##  Change to upper case like DESI
        for colname in zbest.colnames:
          if colname.islower():
            zbest.rename_column(colname, colname.upper())

        ##  Add brickname column
        zbest['BRICKNAME'] = meta['BRICKNAME']

        write_zbest(os.environ['LBGCMB'] + '/redrock/dat/zbest/%s/%s/%s.fits' % (survey, subpath, output), zbest, overwrite=True)
    
        run_time = time.time() - start_time
    
        json.dump(meta, open(os.environ['LBGCMB'] + '/redrock/dat/meta/%s/%s/%s_meta.json' % (survey, subpath, output), 'w'))    

        print('\n\n\nTime: %.2lf s.  Percentage complete: %.2lf;  Redshift: %.2lf, Exposure: %.2lf \n\n\n' % (time.time() - start_time, 100. * (jj + ii * len(redshifts)) / nruns, redshift, exposure))
        
    print("\n\nDone.\n\n")
