import  os
import  numpy              as      np
import  desisim
import  desisim.templates
import  matplotlib.pyplot  as      plt
import  pylab              as      pl  
import  astropy.units      as      u 
import  astropy.constants  as      const

from    desitarget.cuts    import  isELG_colors, isLRG_colors, isBGS_colors, isQSO_colors
from    scipy.signal       import  medfilt
from    astropy.io         import  fits
from    astropy.table      import  Table, Column, hstack
from    numpy.random       import  uniform


ngal               =     250
dband              =      'g'
pure_noise         =       0
target_type        =    'elg'
mag_degrade        =    True
mag_drop           =     4.0
plot               =   False
save               =    True

metaname           =    os.environ['BEAST'] + '/gal_maker/dat/Tables/galmaker-%s-meta.txt' % target_type
infile             =  os.environ['SCRATCH'] + '/desi/simspec/%s-input-spectra.fits' % target_type


print('\n\nWelcome to gal. maker.\n\n')

target_types       =  {'qso': desisim.templates.QSO, 'lrg': desisim.templates.LRG, 'elg': desisim.templates.ELG,\
                       'bgs': desisim.templates.BGS}

##  rlfux cut in desitarget.cuts.  Flux in nanomaggies.  22.39 mag.
cut_types          =  {'qso': isQSO_colors,          'lrg': isLRG_colors,          'elg': isELG_colors,\
                       'bgs': isBGS_colors}                                                                                     

##  https://github.com/desihub/tutorials/blob/master/simulating-desi-spectra.ipynb;  ## (rflux=0.07)
gal_maker          =  target_types[target_type](normfilter_north='BASS-' + dband, normfilter_south='decam2014-' + dband,\
                                   colorcuts_function=(cut_types[target_type]))

if target_type == 'elg':
  ##  2D array of flux [1e-17 erg/s/cm2/A]  
  ##  flux, wave, meta, objmeta         =  gal_maker.make_templates(nmodel=ngal, redshift=uniform(0.4, 0.6, ngal), mag=uniform(18.0, 19.0, ngal),\
  ##                                                                minoiiflux=8.e-17 * u.erg / u.cm / u.cm / u.s, nocontinnum=True)
  ##  flux, wave, meta, objmeta         =  gal_maker.make_templates(nmodel=ngal, minoiiflux=8.e-17, nocontinuum=True, nocolorcuts=False) 
  flux, wave, meta, objmeta             =  gal_maker.make_templates(nmodel=ngal, minoiiflux=0.0, nocontinuum=False, nocolorcuts=False, magrange=(18.0, 22.0))

else:
  ##  Mostly QSO.
  flux, wave, meta, objmeta = gal_maker.make_templates(nmodel=ngal, magrange=(17.0, 23.4))

if mag_degrade:
  ##  Drop mags. and fluxes by 2 mags.  Remove other columns.
  ##  Note:  ignores correlation between e.g. doublet flux and continuum.
  ##         leave equivalent widths alone.
  ##  objmeta['OIIFLUX']    *= 10. ** (-2. / 2.5)
  ##  objmeta['HBETAFLUX']  *= 10. ** (-2. / 2.5)
  meta['MAG']  -=  mag_drop
  ffactor       =  10. ** (-mag_drop / 2.5)

  flux         *=  ffactor
  
  for x in ['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2']: 
    meta[x]    *=  ffactor

  for x in ['OIIFLUX', 'HBETAFLUX']:
    objmeta[x] *=  ffactor

  for x in ['D4000', 'OIIDOUBLET', 'VDISP', 'OIIIHBETA', 'OIIHBETA', 'NIIHBETA', 'SIIHBETA']:
    objmeta.remove_column(x)

print(meta)
print(objmeta)

if target_type != 'qso':
  ##  Can ignore qso objdata.  As list, PCA coeff throw error.                                                                                                                                                                                                                                
  meta               =   hstack([meta, objmeta], join_type='outer', metadata_conflicts='warn')

if pure_noise:
  target_type        =  'purenoise'
  flux               =   np.zeros_like(flux)

else:
  meta['g']          =  Column(data = 22.5 - 2.5 * np.log10(np.array(meta['FLUX_G'])), name='g')  ##  Astropy Table column.
  meta['r']          =  Column(data = 22.5 - 2.5 * np.log10(np.array(meta['FLUX_R'])), name='r')  ##  Astropy Table column. 
  meta['z']          =  Column(data = 22.5 - 2.5 * np.log10(np.array(meta['FLUX_Z'])), name='z')  ##  Astropy Table column.

  for band in ['g', 'r', 'z']:
    print('\nMag limits:  %.4lf < %s < %.4lf' % (np.min(np.array(meta[band].quantity)[np.isfinite(np.array(meta[band].quantity))]), band, np.max(np.array(meta[band].quantity)[np.isfinite(np.array(meta[band].quantity))])))

## print('wave',                   wave)
## print('flux.shape',       flux.shape)

## print('\n\nUnique redshifts: ', meta['REDSHIFT'])
## print('\n\nUnique redshifts: ', meta['g'], meta['r'], meta['z'])

'''
##  -- OII Line flux calc. --  ##
dwave      =  wave[1] - wave[0]

for i, xx in enumerate(flux): 
  pl.clf()

  redshift     = meta['REDSHIFT'][i]
  OIIa         = 3727.092 * (1. + redshift)
  OIIb         = 3729.875 * (1. + redshift)
  OII          = 3728.483 * (1. + redshift)
  dlambda      =      5.1 * (1. + redshift)

  linecut      = ((OII - dlambda) <= wave) & (wave <= (OII + dlambda))
  continuumcut = ((OII + dlambda) <= wave) & (wave <= (OII + dlambda + 1.))

  wavecut      = wave[linecut]
  fluxcut      =   xx[linecut] 
  fluxcut     *= 1.e-17
  lineflux     = np.sum(fluxcut) * dwave
  contflux     = np.sum(xx[continuumcut]) * dwave
  contflux     = np.sum(np.ones_like(fluxcut)) * np.mean([fluxcut[0], fluxcut[-1]]) * dwave

  lineflux    -= contflux

  print('%.6le' % (np.sum(fluxcut) * dwave))
  
  pl.axvline(x=OIIa,     ymin=0., ymax=1., c='k')
  pl.axvline(x=OIIb,     ymin=0., ymax=1., c='k')
  pl.axhline(y=np.mean([fluxcut[0], fluxcut[-1]]), xmin=0., xmax=1., c='k')

  pl.plot(wavecut, fluxcut, 'o', markersize=4, label=str(i))
  pl.title('\t%.2le, %.2le [ergs/s/cm/cm] \t %.1lf per cent' % (lineflux, np.float(meta['OIIFLUX'][i]), 100. * (lineflux - meta['OIIFLUX'][i]) / meta['OIIFLUX'][i]))

  pl.legend()
  pl.show()
  
  break
'''

if plot:
  ##  --  Sanity check plots --
  pl.clf()

  plt.hist(meta['REDSHIFT'], 20, (0,5))
  plt.xlabel('redshift')

  pl.show()
  pl.clf()
  
  plt.hist(meta['g'] , 20, (15, 25))
  plt.xlabel('g magnitude')

  pl.show()
  pl.clf()
  
  for i, x in enumerate(flux):
    pl.semilogy(wave, x, label=str(meta['r'][i]))

  pl.legend(ncol=4)
    
  pl.show()

if save:
  meta.write(metaname, format='ascii', overwrite=True)

  hdr            = fits.Header()

  hdr['EXTNAME'] = 'WAVELENGTH'
  hdr['BUNIT']   = 'Angstrom'

  fits.writeto(infile, wave, header=hdr, overwrite=True)

  hdr['EXTNAME'] = 'FLUX'
  hdr['BUNIT']   = '10^-17 erg/(s*cm^2*Angstrom)'  # Satisifes FITS standard AND Astropy-compatible.

  fits.append(infile, flux, header=hdr)

  print('Written to %s.' % infile)

print('\n\nDone.\n\n')
