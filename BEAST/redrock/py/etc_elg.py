import  os
import  pickle
import  pylab                   as     pl
import  numpy                   as     np
import  matplotlib.pyplot       as     plt
import  matplotlib              as     mpl
import  matplotlib.colors       as     colors

from    astropy.table           import Table, vstack
from    utils                   import latexify
from    redrock.results         import read_zscan
from    mpl_toolkits.axes_grid1 import make_axes_locatable


if __name__ == "__main__":
  print("\n\nWelcome to exposure.\n\n")

  dband           =     'g'
  target          =   'elg'
  survey          =  'desi'
  base_exp        =   1500.        ##  [Seconds].
  repeat          =      1
  nexposures      =     10         ##  Scaling,  not coadd.
  catch_outliers  =  False
  catch_tempdegen =  False

  printit         =   True
  ddir            =     ''

  root            = os.environ['BEAST']

  ##  Truth values.
  names           = ['TARGETID_1', 'OBJTYPE', 'SUBTYPE', 'TEMPLATEID', 'SEED', 'REDSHIFT', 'MAG', 'MAGFILTER',\
                     'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'TARGETID_2',\
                     'OIIFLUX', 'HBETAFLUX', 'EWOII', 'EWHBETA', 'D4000', 'VDISP', 'OIIDOUBLET', 'OIIIHBETA', 'OIIHBETA', 'NIIHBETA', 'SIIHBETA', 'g', 'r', 'z']

  Truth           = os.environ['BEAST'] + '/gal_maker/dat/Tables/galmaker-elg-meta.txt'
  Truth           = Table.read(Truth, format='ascii', names=names)
  Truth           = Table([Truth['OBJTYPE'], Truth['SUBTYPE'], Truth['REDSHIFT'], Truth['EWOII'], Truth[dband]],\
                          names=('OBJTYPE', 'SUBTYPE', 'REDSHIFT', 'EWOII', dband))
  
  Truth           = vstack([Truth] * repeat)

  ##  No redshift success for any of the exposures.
  exptimes        =  np.ones_like(Truth['REDSHIFT'].quantity) * 1.e99
  
  ##  Digitize in EW and replace by bin centre.
  bins            =  np.linspace(np.array(Truth['EWOII'].quantity).min() - 1., np.array(Truth['EWOII'].quantity).max() + 1., 4)
  inds            =  np.digitize(np.array(Truth['EWOII'].quantity), bins)
  cens            =  0.5 * (bins + np.roll(bins, 1))
  cens[0]         =  np.NaN
  cens            =  np.round(cens, decimals=2)
  
  Truth['EWOII']  =  cens[inds]

  TRUEZS          =  np.array(Truth['REDSHIFT'].quantity).astype(np.float)
  TRUEEWS         =  np.array(Truth['EWOII'].quantity).astype(str)
  TRUEMAGS        =  np.array(Truth[dband].quantity).astype(np.float)

  EWs             = [str(x) for x in cens]
  results         =  dict(zip(EWs, [[] for x in EWs]))

  for scaling in np.arange(1, nexposures, 1):
    exposure       =  scaling * base_exp            ##  Seconds.
    
    print('\nSolving for exposure: %d' % exposure)

    infile         =  os.environ['CSCRATCH'] + '/desi/simspec/%s/elg-desi-spectra-exp-%d.fits'    % (ddir, exposure)  
    output         =  infile.split('/')[-1].split('.fits')[0]

    rrh5file       =  os.environ['CSCRATCH'] + '/desi/simspec/%s/%s-%s-spectra-exp-%d-rr.h5'      % (ddir, target, survey, exposure)
    zbestfile      =  os.environ['CSCRATCH'] + '/desi/simspec/%s/%s-%s-spectra-exp-%d-zbest.fits' % (ddir, target, survey, exposure)
    
    Data           =  Table.read(zbestfile)
      
    ##  print
    ##  print(Data)
    ##  print

    (zscan, zfit)  =  read_zscan(rrh5file)

    ##  Unique full spectral types ('spectype' + 'sub_type')                                                                   
    utypes = set([x + '-' + zfit['subtype'][kk] if zfit['subtype'][kk] != '' else zfit['spectype'][kk] for kk, x in enumerate(zfit['spectype'])])

    ##  zfit.info
    ##
    ##  set(zfit[:]['spectype'])
    ##  {'GALAXY', 'LBG', 'QSO', 'STAR'}
    ##
    ##  set(zfit[:]['subtype'])
    ##  {'', 'A', 'B', 'F', 'K', 'M'}
    
    zz         =  zfit[zfit['znum'] == 0]          ##  Results for best-fitting redshift of each target.                                                    

    ids        =  zz['targetid']
    zbests     =  zz['z']
    zerrs      =  zz['zerr']
    minchi2s   =  zz['chi2']                       ##  Min. chi2 calculated at finer z resolution for expected target type. 
                                                   ##  Previously, zscan[kk][template_type]['zchi2'].min().    
    zwarns     =  zz['zwarn']

    ##  Run through targets and set success to 0 if redrock warning or true z and fitted z differ by more than error.  
    for kk, id in enumerate(ids):        
      ##  Rank of best fitting redshift.
      zbest      =  zbests[kk]
      zerr       =  zerrs[kk]
      minchi2    =  minchi2s[kk]                   ##  Chi2 calculated at finer z resolution for expected target type.  
                                                   ##  Previously, zscan[kk][template_type]['zchi2'].min().
      zwarn      =  zwarns[kk]
      targetid   =  ids[kk]

      ##  True z.
      truez      =  TRUEZS[kk]
      trueEW     =  TRUEEWS[kk]      
      truemag    =  TRUEMAGS[kk]

      ##  Digitize equivalent width.
      trueEW     =  np.round(cens[np.digitize(trueEW, bins)], decimals=2)

      ##  For each target, find if best-fitting redshift had a warning from redrock, or if the fitted redshift is 5. * zerr from the truth. 
      if (zwarn != 0):
        success  = 0
    
        if printit:
          print('\nTarget ID %d caught by ZWARN: z= %.3lf, m=%.3lf, exp=%d; z best = %.3lf +- %.3le, X2 = %.3lf, warning: %s, success: %d' % (targetid, truez, truemag,\
                                                                                                                                              exposure, zbest,\
                                                                                                                                              zerr, minchi2,\
                                                                                                                                              zwarn, success))
          
      ##  Catch catastrophic outliers.
      elif catch_outliers & (zwarn == 0) & (np.abs(truez - zbest) > 1. * zerr):
        success   = 0

        if printit:
          print('\nTarget ID %d caught by ZOUTLIER: z= %lf, m=%.3lf, exp=%d; z best = %lf +- %.6le, X2 = %.3lf, warning: %d, success: %d' % (targetid, truez, truemag,\
                                                                                                                                             exposure, zbest,\
                                                                                                                                             zerr, minchi2,\
                                                                                                                                             zwarn, success))
      
      if catch_tempdegen & (success == 1):
        ##  Now check on successful discrimination from other templates.   
        ##  Already implemented by redrock.
        ##
        for spectype in zscan[kk].keys():
          if spectype != 'LBG':
            zx      = zscan[kk][spectype]
            rchi2   = (zx['zchi2'] + zx['zchi2']).min()
            rchi2_z = zx['redshifts'][(zx['zchi2'] + zx['zchi2']) == rchi2]

            ##  Difference in chi^2 must be at least twenty.
            ##  Note:  Penalty?
            if (rchi2 - minchi2) > 20:
              continue

            else:
              success = 0
            
              print('\nTarget ID %d caught by TEMPLATE DEGENERACY: z= %lf, m=%.3lf, exp=%d; z best = %lf +- %.6le, X2 = %.3lf, confusion: %s at z=%.3lf with X2: %.3lf'\
                    % (targetid, truez, truemag, exposure, zbest, zerr, minchi2, spectype, rchi2_z, rchi2))
      
      else:
        success   = 1

      if success == 1:
        ##  Ok, good redshift.  Save to the list.
        print('Adding target ID: %d, z=%.3lf, r=%.3lf, exposure= %.2lf to successful exposures.' % (targetid, truez, truemag, exposure))

        exptimes[kk] = np.array([exptimes[kk], exposure]).min()

  pickle.dump(exptimes, open(root + '/redrock/pickle/%s/%s/exptimes.pkl' % (survey, target), 'wb'))    
  
  print("\n\nDone.\n\n")
