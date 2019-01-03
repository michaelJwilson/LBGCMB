import  pylab                                as      pl
import  numpy                                as      np
import  matplotlib                           as      mpl
import  matplotlib.pyplot                    as      plt
import  astropy.units                        as      u
import  astropy.constants                    as      const

from    astropy.convolution                  import  convolve, Box1DKernel
from    stellar_spectra.stellar_spectra      import  get_onine
from    read_ised                            import  read_ised
from    em_lines                             import  plot_emlines
from    prep_filters                         import  prep_filters
from    redshift_spectra                     import  redshift_spectra
from    utils                                import  latexify


latexify(fig_width=None, fig_height=None, columns=2, equal=False)

def load_composite(mJy=False, getcontinuum=False, npix=64, wave_ext=0.):
    '''                                                                                                                                                  
    Plot the composite LBG spectra from Shapley (2003);                                                                                                   
    https://arxiv.org/pdf/astro-ph/0301230.pdf                                                                                                                    
    See the 2016 update by Steidel ++                                                                                                                            https://arxiv.org/pdf/1605.07186.pdf                                                                                                                         '''

    data         =  np.loadtxt('../quickspectra/spectra/spec-shapley-z0.000.dat')

    mlambda      =  data[0,0]

    dlambda      =  data[1,0] - data[0,0]
    nlambda      =  len(data[:,0])

    ls           =  data[:,0]
    vs           = (ls * u.AA).to(u.Hz, equivalencies = u.spectral())

    Fv           =  data[:,1]
    
    if mJy:
      Fv        *= 1e29  ## [micro Jy].

    if getcontinuum:
        npix     = npix

        ##  Reflection at the far edge for convolution; zero padding on blue end is probably fine.                                                     
        ls       = np.concatenate([ls, ls[-1] + ls])
        Fv       = np.concatenate([Fv, Fv[::-1]])

        Fv       = convolve(Fv, Box1DKernel(npix))

        ## Cut reflection range. 
        ls       = ls[:nlambda]
        Fv       = Fv[:nlambda]

        vs       = (ls * u.AA).to(u.Hz, equivalencies = u.spectral())

    if wave_ext > 0.:
        ## Blue extension to 0.0 angstroms and red extension to wave_ext.
        ext      = np.arange(0.0, mlambda, dlambda)

        ls       = np.concatenate([ext, ls])
        Fv       = np.concatenate([np.zeros_like(ext), Fv])

        nir_amp  = Fv[(ls > 1600.)].mean()
        nir_std  = Fv[(ls > 1600.)].std()

        ext      = np.arange(2000. + dlambda, wave_ext + dlambda, dlambda)

        ls       = np.concatenate([ls,   ext])
        Fv       = np.concatenate([Fv, nir_amp * np.ones_like(ext) + nir_std * np.random.uniform(low=-.5, high=0.5, size=len(ext))])

        vs       = (ls * u.AA).to(u.Hz, equivalencies = u.spectral())

    return  vs.value, ls, Fv


if __name__ == "__main__":        
    print('\n\nWelcome to composite.\n\n')

    Age                     =  25.
    redshift                =  0.001

    ## Age [Myr]; vs [Hz]; SEDS [ergs/s/Hz]; ls [A]; Ll [ergs/s/AA].
    ages, vs, Fv, ls, Ll    =  read_ised('GALAXEV/models/Padova1994/salpeter/bc2003_hr_m72_salp_ssp.ised', AgeMyr = Age, printit = False)

    vs, Fv                  =  redshift_spectra(vs, Fv, redshift=redshift)
    ls                      = (1.e10 / vs) * const.c.to('m/s').value                     ##  Redshifted wavelengths.

    pl.plot(ls, Fv / Fv.max(), c='k', alpha=1.0, label=r'$z$=%.1lf' % redshift + ';  %.1lf Myr' % ages)
    
    ## Shapley rest-frame composite.
    vs, ls, Fv              =  load_composite(mJy=False, getcontinuum=False, npix=64, wave_ext=0.)

    vs, Fv                  =  redshift_spectra(vs, Fv, redshift=redshift)
    ls                      = (1.e10 / vs) * const.c.to('m/s').value                     ##  Redshifted wavelengths. 
    
    pl.plot(ls, 1.4 * Fv / Fv.max(), 'g', label=r'14.3$\AA$ composite', alpha=0.85)
    
    ## ls, Fv               = get_onine(mJy = False)                                     ##  Spectra of an O9 star.
    ## pl.plot(ls, spectra * factor, 'r')
    
    ##  LSST filters.
    filters                 = prep_filters(filter = 'LSST')
        
    for i, band in enumerate(filters.keys()):
        if filters[band]['ppkey'][-2] is not 'y':          
          pl.fill_between(filters[band]['ls'], filters[band]['Ts'] * 0.85, alpha=0.9, label=filters[band]['ppkey'])

    ##  VIDEO filters.
    filters                 = prep_filters(filter = 'VIDEO')

    for i, band in enumerate(filters.keys()):
        if band not in ['H', 'K']:
          pl.fill_between(filters[band]['ls'], filters[band]['Ts'] * 0.40, alpha=0.9, label=filters[band]['ppkey'])
    
    pl.xlim(0.0, 1.2e4)
    pl.ylim(0.0,   1.2)
    
    ## plot_emlines(ymax=0.6)
    
    pl.xlabel(r"$\lambda \ [\AA]$")
    
    ## pl.ylabel(r"$F_{\nu} \ [\mu \rm{Jy}]$")
    pl.ylabel(r"$F_{\nu}$ (arbitrary normalisation)")     

    pl.legend(ncol=2, loc=1)
    
    ## Limited to 2000 A.
    pl.savefig("plots/composite.pdf", bbox_inches='tight')
