import  numpy              as      np
import  matplotlib.pyplot  as      plt

from    prep_Llls          import  prep_Llls
from    pmh                import  Pmm, get_PkInterps
from    prep_camb          import  CAMB
from    utils              import  comoving_distance
from    collections        import  OrderedDict
from    pz2nbar            import  nbar_convert
from    whitebook_pz       import  whitebook_pz, Ng, lsst_chang_pz
from    utils              import  latexify, sci_notation
from    scipy.interpolate  import  interp1d
from    schmittfull_bz     import  desi_eff
from    growth_rate        import  growth_factor
from    pz_tools           import  percentiles
from    Cgg                import  Cgg


##  -- TO DO. -- 
## 

##    eqn. (4) of D1.1, pg. 46 of latest LSST WB, with kmax = 0.3
KMAX      = 0.3

def get_desiCz():
    ##  Saved file of zmax, N(z < zmax) [gals. per sq. deg.], eff. b(zmax).
    ##  Produced by desi_pz.py
    data  = np.loadtxt('dat/effective_desi.dat')
    
    zs    = data[:,0]
    cs    = data[:,1]

    nbar  = cs[-1]
    cs    = cs / cs[-1]

    pz    = cs - np.roll(cs, 1)
    pz[0] = cs[0]

    ##  Shift by half the bin width. 
    zs   -= (zs[1] - zs[0]) / 2. 

    pz    = interp1d(zs, pz, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

    return  nbar, pz

def lenses_bz(zee):
    aa        = 1. / (1. + zee)

    return  0.95 / growth_factor(aa)
    
def Cij(Pk_interps, Llls, zmin, zmax): 
    zee       = (zmin + zmax) / 2.

    chi       =  comoving_distance(zee)
    dchi      =  comoving_distance(zmax) - comoving_distance(zmin)

    ##  Extended Limber approximation.
    ks        =  (Llls + 0.5) / chi               ## For the Phh evaluation in the integral, we take a zeff approx.      
                                                  ## i.e. \int dz .... Phh(zeff).                                                
    
    ##  Effectively measures the non-linear matter power spectrum, but normalised down (by flux-like factor).   
    return  Pmm(Pk_interps, ks, zee, 'nlinear') / chi ** 2. / dchi
    
def Fisher(Pk_interps, Llls, zs, tNs, tNp, pz, bz, fsky=0.1, fover=0.0, sources = False, printit=True):
    '''    
    Returns the variance on \hat Np = F^{-1}_ii.
    '''

    result              =  OrderedDict()

    for i, zee in enumerate(zs):      
      LMAX              =   np.ceil(KMAX * comoving_distance(zee) - 0.5).astype('int')
      LCUT              =  np.where(Llls <= LMAX)

      ##  Returns linear bias of spectroscopic galaxies at each redshift. 
      bs    =  bz(zee)

      ##  Returns number of spectroscopic redshifts in each shell.
      Ns    =  nbar_convert(tNs, unit='str') * pz(zee) 
      Ns   *=  dz

      if sources:
        ##  Assumed LSST Whitebook for sources, b(z) = 1. + z.                                                                           
        bp  =  1. + zee

      else:
        ## Lenses, pg. 47 of latest whitebook.                                                                                           
        bp =   lenses_bz(zee)

      ##  Returns number of photometric redshifts in each shell.
      Np   =   nbar_convert(tNp, unit='str') * lsst_chang_pz(zee, ilim=25.3, source = False)
      Np  *=   dz
      
      ## Returns shot noise: wp = Np for fsat = 0.0; eqn. (11) of MQW16. 
      wp   =   Np
      ws   =   Ns

      wps  =   fover * np.min([wp, ws])

      cgg  =   Cij(Pk_interps, Llls, zee - dz / 2., zee + dz / 2.)
      cgg  =   Cgg(Pk_interps, Llls, zee - dz / 2., zee + dz / 2., pz, bz,\
                   survey_pz2 = lambda z: lsst_chang_pz(z, ilim=25.3, source=sources), bz2=lenses_bz, zeff=True)   

      ## eqn. (44); sum over k has implicit sum over ell and m.                                                                                 
      result[zee]           = {'Ns': Ns, 'Np': Np, 'bp': bp, 'bs': bs, 'cgg': cgg, 'wp': wp, 'num': (bp * bs * cgg) ** 2.}

      result[zee]['A00z']   = cgg * (bp * Np)**2. + wp        ## eqn. (22),  < p   p   >_L
      result[zee]['A0i']    =  bp * Np * bs * Ns * cgg + wps  ## eqn. (23),  < p   s_i >_L
      result[zee]['Aii']    =    (bs * Ns) ** 2. * cgg + ws   ## eqn. (24),  < s_i s_j >_L  
      result[zee]['A0i_j']  =       bp * bs * Ns * cgg        ## eqn. (25),  < p   s_i >_L, j 

      ## eqn. (43); to be later normalised. 
      result[zee]['beta']   = (Np * bp) ** 2. * cgg 
    
    ## Total spec. zs over the whole redshift range. 
    nspec   = 0.0

    ## ... and normalisation over redshift for beta.                                                                                   
    bnorm   = np.zeros_like(result[zee]['beta'])
    A00     = np.zeros_like(result[zee]['A00z'])   ## < p * p >_L
    sshift  = np.zeros_like(result[zee]['A0i'])

    for zee in result:
        nspec  += result[zee]['Ns']
        
        ## Array for each L value.
        A00    += result[zee]['A00z'] 
        bnorm  += result[zee]['beta']
        
        ## Sum to S, eqn. (29);  Together with A00, yields S.
        sshift += result[zee]['A0i'] ** 2. / result[zee]['Aii'] 
    
    for zee in result: 
        ## Schur-Limber limit, eqn. (44).
        result[zee]['beta']         /=  bnorm

        ## Fractional error from eqn. (44);  Schur-Limber limit of small r (due to shotnoise, or redshift overlap). 
        result[zee]['ratio']         =  result[zee]['num'] / A00

        ## Diagonal by definition in this limit;  Impose LMAX cut.  
        result[zee]['Fii']           =  result[zee]['Ns'] * fsky * np.sum((2. * Llls[LCUT] + 1.) * result[zee]['ratio'][LCUT])
        result[zee]['var_ii']        =    1. / result[zee]['Fii']

        ## NOTE:  the error is for b * N(z), this fractional error is screwey unless renormalised almost immediately.
        result[zee]['ferr_ii']       = np.sqrt(result[zee]['var_ii']) / result[zee]['Np']
    
    ## Limber approximation, but not in the Schur limit (small cross-correlation).    
    for zi in result:
        ## Schur(L)
        result[zi]['S']              = A00  / (A00 - sshift)
        result[zi]['r']              = result[zi]['A0i'] / np.sqrt(A00 * result[zi]['Aii'])
        
    ## Non-diagonal Fisher matrix outwith the Schur limit. 
    Fisher                           = np.zeros((len(zs), len(zs)))

    for i, zi in enumerate(result):
      Interim                        =  (2. * Llls + 1) * result[zi]['S'] * result[zi]['A0i_j'] * result[zi]['A0i_j'] / result[zi]['Aii'] / A00
      Fisher[i, i]                   =  fsky * np.sum(Interim[LCUT])

      for j, zj in enumerate(result):
        Interim                     +=   (2. * Llls + 1) * fsky * 2. * result[zi]['S'] ** 2. * result[zi]['r'] * result[zj]['r'] \
                                             * np.sqrt(1. / result[zi]['Aii'] / result[zj]['Aii']) * result[zi]['A0i_j'] * result[zj]['A0i_j'] / A00  

        Fisher[i, j]                +=  np.sum(Interim[LCUT])

    iFish    = np.linalg.inv(Fisher)
    diFish   = np.diag(iFish)

    for i, zz in enumerate(result):
        dstr  = "\tz:  %.2lf \t\t Schur-Limber:  %.2lf \t\t Limber:  %.2lf" % (zz, 100. * result[zz]['ferr_ii'], 100. * np.sqrt(diFish[i])/result[zz]['Np'])
        dstr += "\t\t Schur at Lmin: %.2lf, and Lmax: %.2lf"                % (result[zi]['S'][0], result[zi]['S'][-1])

        if printit:
          print  dstr
    
    ##  Define output. 
    output = []
    
    for outz in zs: 
        index   = np.where(np.abs((zs - outz)) == np.min(np.abs(zs - outz)))[0]

        ## Save 
        output.append([outz, result[outz]['Ns'], result[outz]['Np'], 100. * result[zs[index][0]]['ferr_ii'],\
                                                                     100. * np.sqrt(diFish[index]) / result[zs[index][0]]['Np']])
    
    output = np.array(output)
    
    if printit:
      print  output
    
    return   output


if __name__ == '__main__':
  import  os
  import  pylab              as      pl

  from    pickle             import  dump, load
  from    pz_tools           import  percentiles  


  print('\n\nWelcome to a McQuinn and White clustering redshift forecaster.')
  
  dz           =  0.1    
  zs           =  np.arange(dz, 1.6, dz)

  fover        =  0.00
  fsky         =  5060. / 41253.
    
  Ns, pz       =  get_desiCz()
  Np           =  Ng(25.3, fmask=0.12)

  evaluate     =  True
  sources      =  False

  if evaluate:
    ## Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                    
    cambx        =  CAMB()
    Pk_interps   =  get_PkInterps(cambx)

    ## Input NLlls is ignored in the log10=False case.                                                                                   
    NLlls, Llls, nmodes =  prep_Llls(NLlls = 60, Lmin = 60., Lmax = 5000., log10=False)

    if sources: 
      ##  Determining the source distribution:                                                                                                                     
      pz      = lambda z:  lsst_chang_pz(z, ilim=25.3, source = True)
      slices  = percentiles(pz, printit=False, asarray=False)

      zs      = np.array([0.0] + slices + [10.0])
      result  = Fisher(Pk_interps, Llls, zs, Ns, Np, pz, desi_eff, fsky=fsky, fover=fover, sources=sources)

    else:
      result  = Fisher(Pk_interps, Llls, zs, Ns, Np, pz, desi_eff, fsky=fsky, fover=fover, sources=sources)
      np.savetxt('dat/mqw_cadence.txt', result, fmt='%.6le', delimiter='\t')


  ##  And plot.
  latexify(columns=2, equal=False, fontsize=10, ratio=None, ggplot=True)

  color  = 'dodgerblue'
  result = np.loadtxt('dat/mqw_cadence.txt')
  
  pl.plot(result[:,0], result[:,2], c=color)

  ind = np.where(result[:,4] < 100.)[0]

  pl.errorbar(result[ind,0], result[ind,2], result[ind,2] * result[ind,4] / 1.e2, label='LSST Lenses', c=color)

  pl.fill_between(result[ind,0], result[ind,2] * (1. - result[ind, 4] / 1.e2), result[ind, 2] * (1. + result[ind, 4] / 1.e2), color=color, alpha=0.5)

  pl.xlabel(r'$z$')
  pl.ylabel(r'$N_p(z)$')

  pl.xlim(0.0, 1.7)

  pl.legend(loc=1)
  
  pl.savefig('plots/mqw_cadence.pdf', bbox_inches='tight')
  
  print('\n\nDone.\n\n')
