import  os
import  numpy              as  np
import  pylab              as  pl
import  matplotlib.pyplot  as  plt
import  threedhst.eazyPy   as  eazy
import  threedhst.utils    as  utils


ROOT       = os.environ['SCRATCH']
OUTPUT     = ROOT + '/MUSE/Rafelski15/EAZY'

LBGCOLORS  = ['u-g', 'g-r', 'r-i', 'i-z']
results    = []

for idx in np.arange(9970)[:9960]:
  ##  eazyPy.py, symlinked in current directory. 
  result = eazy.plotExampleSED(idx=idx, pdfname='plot.pdf', MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY=OUTPUT,\
                               CACHE_FILE='Same', lrange=[1.e3, 1.e5], axes=None, individual_templates=False, fnu=False)
  
  results.append(result)
            
  if idx % 100 == 0:
    np.savetxt('colors_mag_25p0_1k.txt', np.array(results), fmt='%.6lf')

np.savetxt('colors_mag_25p0_1k.txt', np.array(results), fmt='%.6lf')              

'''
################################################
#### Pull out data from the BINARY_OUTPUTS files
################################################

## SED, data & fit
sed = eazy.getEazySED(idx=17, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY=OUTPUT, CACHE_FILE='Same',\
                      scale_flambda=True, verbose=False, individual_templates=False)

lambdaz, temp_sed, lci, obs_sed, fobs, efobs = sed

axes[0].scatter(lci, obs_sed, color='orange', zorder=2)
axes[0].scatter(lci,    fobs, color='green',  marker='s', s=150, zorder=2)

## p(z)
zgrid, pzi, prior = eazy.getEazyPz(17, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY=OUTPUT,
                                   CACHE_FILE='Same', binaries=None, get_prior=True)
                                   
axes[1].plot(zgrid, pzi, color='black')
axes[1].plot(zgrid, prior/prior.max()*pzi.max(), color='purple')

plt.savefig('eazy_fit_2.png')


################################################
#### Investigate template fit residuals for zeropoints
################################################
eazy.show_fit_residuals(root='photz', PATH=OUTPUT, savefig='fit_residuals.png', adjust_zeropoints='zphot.zeropoint',\
                        fix_filter=None, ref_filter=28, get_resid=False, wclip=[1200, 30000.0])

## Run again using the zphot.zeropoint file you just made, this can be done iteratively
params                   = eazy.EazyParam('zphot.param.default')
params['GET_ZP_OFFSETS'] = 'y'

params.write('zphot.param')

os.system('../src/eazy -p zphot.param -z zphot.zeropoint')

## Check z_phot vs z_spec again
## eazy.zPhot_zSpec('OUTPUT/photz.zout')
plt.savefig('zphot_zspec_2.png')
'''

print("\n\nDone.\n\n")
