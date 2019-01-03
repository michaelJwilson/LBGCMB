import  numpy              as      np
import  healpy             as      hp
import  matplotlib.pyplot  as      plt

from    healpy.pixelfunc   import  npix2nside, nside2npix
from    astropy.table      import  Table


def  load_hp(fpath='cmb/dat/mask/mask_40pc.fits', view=False, printit=False):
  extension = fpath.split('.')[-1]

  if extension == 'fits':
    from    astropy.io  import  fits

    ##   Survey files of the form:  https://lambda.gsfc.nasa.gov/toolbox/footprint/configfile.cfm
    T      = Table.read(fpath)

    ORDER  = T.meta['ORDERING']      
    NROWS  = len(T)      

    NSIDE  = T.meta['NSIDE']
    NPIX   = nside2npix(NSIDE)

    count  = 0.0

    mask   = []

    for i in np.arange(NROWS):
      mask  += [element for element in T[i][0]]

    mask   = np.array(mask)

  else:
    mask   =  np.loadtxt(fpath)

    NPIX   =  len(mask)                                                                                                                                                                               
    NSIDE  =  npix2nside(NPIX)
  
  print('Loaded healpy mask: %s as %s (NPIX:  %d, \t NSIDE:  %d)' % (fpath.split('/')[-1], extension, NPIX, NSIDE))

  if printit:
    print('Min:  %s;  Max:  %s' % (mask.min(), mask.max()))
  
  ##  Element check.
  mask[mask > 1.0] = 1.0
  mask[mask < 0.0] = 0.0
  
  return  mask

def print_grid(surveys, result):
  print('\n\n\n\n')

  print(''.join('\t %s'.ljust(10) % survey for survey in surveys))

  for kk, row in enumerate(result):
    print(''.join('%s' % surveys[kk]) + ''.join('\t %.2lf'.ljust(10) % crossarea for crossarea in row))

def plots(maps, surveys, clf=True, save=True, title=None, pure_noise = True):
  import  pylab             as      pl

  from    utils             import  latexify
  

  latexify(columns=2, equal=False, fontsize=10, ratio=None, ggplot=True)

  for kk, survey in enumerate(surveys):
    if survey in ['Planck', 'SO', 'AdvACT']:
      fmap  = 1.e6 * cmb_map(survey, pure_noise = pure_noise)

    else:
      fmap  = np.ones_like(maps[kk]) 

    if clf:
      pl.clf()

    ##  C equates to equatorial. ## min=-150., max=150.
    hp.visufunc.mollview(maps[kk] * fmap, coord=['C', 'C'], flip='astro', rot=(0., 0., 0.), hold=False, format='%.1lf', cbar=True,\
                         unit=r'$\mu K$', cmap='RdBu_r')
    
    hp.graticule(dpar = 15., dmer = 30.)

    if save:
      if title is None:
        pl.title('%s in Mollweide'          % survey)
        pl.savefig('plots/footprint/%s.pdf' % survey)

      else:
        pl.title('%s' % title)
        pl.savefig('plots/footprint/%s.pdf' % ('_' + title))

def rotate_map(map, inco='C', outco='G'):
    nside   = hp.npix2nside(len(map))
    
    ##  Get theta, phi for non-rotated map.
    t, p    = hp.pix2ang(nside, np.arange(len(map)))
    
    ##  Define a rotator
    r       = hp.Rotator(coord=[inco, outco])
    
    ##  Get theta, phi under rotated co-ordinates
    tt, pp  = r(t, p)
        
    return  hp.get_interp_val(map, tt, pp)

def cmb_map(cmbexp = 'Planck', nside=512, pure_noise = True, plot_Cls = False, plot_map = False):
  import  pylab              as      pl 

  from    bolometers         import  bolometers
  from    prep_camb          import  CAMB, Clxy
  from    DetectorNoise      import  DetectorNoise  
  from    prep_Llls          import  prep_Llls
  from    utils              import  prefactor


  fsky, thetab, DeltaT, iterative    =   bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                         bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

  NLlls, Llls, nmodes                =   prep_Llls(NLlls = 60, Lmin = 1., Lmax = 5000., log10=False)
  Cls                                =   DetectorNoise(Llls, thetab, DeltaT, type='TT')

  if not pure_noise:
    ## Prepare pycamb module; linear, non-linear matter P(k) and Cls.
    cambx                              =   CAMB()

    ## No Detector noise -- this is handled by Clxy of prep_camb.                                                                                    
    (lensCl_interps, nolensCl_interps) =   cambx.get_Cls()

    Cls                               +=   lensCl_interps['TT'](Llls)

  map                                  =   hp.sphtfunc.synfast(Cls, nside, lmax=5000, pol=False, pixwin=True, fwhm=0.0, sigma=None, new=False, verbose=True)

  if plot_Cls:
    pl.clf()
    pl.plot(Llls, prefactor(Llls, n=2) * Cls)

    pl.ylim(0.0, 1.e-9)
    
    pl.savefig('Cls.pdf')

  if plot_map:
    pl.clf()

    hp.mollview(map, title='')
    
    pl.savefig('map.pdf')

  return  map

def prep_masks(surveys):
  import  os

  root       =  os.environ['LBGCMB'] + '/mask/'
  
  fpaths     =  {'Planck': 'Planck_MaskInt_UT78_256.fits', 'SPT': 'SPT_150_hits_hpx.fits', 'AdvACT': 'ACT_148_equ_hits_hpx.fits', 'SO': 'mask_40pc.fits',\
                    'QSO': 'BOSS_dr12_qso.fits', 'LSST': 'opsim_nvisits_g.fits', 'HSC': 'hsc_256.txt', 'AdvACT-SOU': 'ACT_148_south_hits_hpx.fits'}

  strategy   =   'deep'

  fpaths['DESI-Y1'] =  'desi/%s/yr0_64.fits' % strategy 
  fpaths['DESI-Y2'] =  'desi/%s/yr1_64.fits' % strategy
  fpaths['DESI-Y3'] =  'desi/%s/yr2_64.fits' % strategy
  fpaths['DESI-Y4'] =  'desi/%s/yr3_64.fits' % strategy
  fpaths['DESI-Y5'] =  'desi/%s/yr4_64.fits' % strategy

  masks             =  []

  for survey in surveys:
      fpath    = root + fpaths[survey]
      masks.append(load_hp(fpath, view=False))

  ## DESI (desi_256.dat) defaults to nested. Convert to ringed.                                                                                                                               
  ## masks[2]  = hp.pixelfunc.reorder(masks[2], n2r=True)                                                                                                                                     

  if 'SO' in surveys:
    ## Convert SO from galactic to ecliptic.                                                                                                                                                  
    masks[surveys.index('SO')]      = rotate_map(masks[surveys.index('SO')], 'C', 'G')

  if 'Planck' in surveys:
    ## Convert Planck from galactic to ecliptic.                                                                                                                                              
    masks[surveys.index('Planck')]  = rotate_map(masks[surveys.index('Planck')], 'C', 'G')

  ## Highest resolution available.                                                                                                                                                            
  nside                             =  np.array([npix2nside(len(mask)) for mask in masks]).max()

  ##  Upgrade all masks to the nside with greatest resoltuion.                                                                                                                                
  masks                             = [hp.pixelfunc.ud_grade(mask, nside) for mask in masks]

  return  masks, nside


if __name__ == '__main__':
    print('\n\nWelcome.\n\n')

    ## 'ACT-SOU', 'SPT', AdvACT
    ## surveys  =  ['Planck', 'SO', 'AdvACT']
    surveys     =  ['Planck', 'SO']
    surveys    +=  ['QSO', 'LSST', 'HSC']
    surveys    +=  ['DESI-Y1', 'DESI-Y3', 'DESI-Y5']

    masks, nside  =  prep_masks(surveys)

    ## Plot survey footprints.                                                                                                        
    ## plots(masks, surveys, clf=False, pure_noise=True, title='')

    ##  Get baseline area for each survey. 
    areas        = [np.sum(mask) * hp.pixelfunc.nside2pixarea(nside, degrees=True) for mask in masks]

    result       = np.array([[np.sum(mask * smask) for smask in masks] for mask in masks])
    result      *= hp.pixelfunc.nside2pixarea(nside, degrees=True)

    print_grid(surveys, result)
    
    print('\n\nDone.\n\n')
