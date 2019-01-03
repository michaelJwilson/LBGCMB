import pickle
import numpy              as      np
import pylab              as      pl
import astropy.constants  as      const

from   numpy.linalg       import  inv
from   scipy.interpolate  import  interp1d
from   params             import  params
from   astropy.cosmology  import  FlatLambdaCDM
from   utils              import  prefactor, comoving_distance
from   Cgg                import  sliced_pz, Ngg
from   scipy.integrate    import  simps
from   scipy.misc         import  derivative
from   dplus              import  growth_factor
from   qso_bz             import  qso_bz

cosmo            = FlatLambdaCDM(H0 = 100.*params['h_100'], Om0 = params['om_m'], Ob0 = params['om_b'])

wiggle           = np.loadtxt('wiggle.dat')
nowiggle         = np.loadtxt('nowiggle.dat')

diff             = wiggle - nowiggle
diff[:,0]        = wiggle[:,0]

diff_interp      = interp1d(    diff[:,0],     diff[:,1], bounds_error=False, fill_value=0.0)
wiggle_interp    = interp1d(  wiggle[:,0],   wiggle[:,1], bounds_error=False, fill_value=0.0)
nowiggle_interp  = interp1d(nowiggle[:,0], nowiggle[:,1], bounds_error=False, fill_value=0.0) 

survey_dzs       = {"LSST": 2.0, "SDSS9": 2.0, "CMASS": 0.15, "QSO": 0.15}  ## Assumed width in dz. 

def lsst_bz(z):
    return 1.0 + z

def qso_bz(z):
    return 0.278 * ((1. + z)**2. - 6.565) + 2.393

def cmass_bz(z):
    return 2.0

def bzs(z, survey):
    if survey == "CMASS":
        return cmass_bz(z)

    elif survey == "QSO":
        return qso_bz(z)

    elif survey == "LSST":
        return lsst_bz(z)

    elif survey == "SDSS9":
        return lsst_bz(z)

    else:
        raise ValueError("Erroneous input to bzs:  ", z, survey)

def Pab(k, ba, bb, alpha = 1.0, sigma = 2., nowiggle=False):
    if nowiggle == False:
        interim = ba * bb * (nowiggle_interp(k / alpha) + diff_interp(k / alpha) * np.exp( - k * k * sigma * sigma / 2.))  ## Include wiggles.
        
        return interim * alpha**2.

    elif nowiggle == True:
        interim = ba * bb *  nowiggle_interp(k / alpha)
        
        return interim * alpha**2.

    else:
        raise ValueError("Erroneous input to Pab.")

def nPab(k, ba, bb, alpha = 1.0, sigma = 2., Aps = 1.0, Am3=0.0, Am2=0.0, Am1=0.0, Ap0=0.0, Ap2=0.0, nowiggle=False):
    ## Pps with additional nuisance parameters to mitigate information emanating from the broad-band power shape.
    result   = Pab(k, ba, bb, alpha, sigma, nowiggle=nowiggle)
    result  *= Aps
    
    result  += Am3 * k **-3.
    result  += Am2 * k **-2.
    result  += Am1 * k **-1.
    result  += Ap0 * k **+0.
    result  += Ap2 * k **+2.

    return result

def Cab(Llls, spec_z, surveya = "LSST", surveyb = "QSO", alpha = 1.0, sigma = 2., zeff=True, nowiggle=False):
  '''
  Angular correlation function: i.e. C_sp(ell, z).
  '''

  z           = np.arange(0.01, 5.0, 0.01)
  chis        = comoving_distance(z)

  ba          = bzs(spec_z, surveya)
  bb          = bzs(spec_z, surveyb)
  
  if zeff == True:
    '''                                                                                                                                                    
    Compute Cgg in the z_eff and first-order Limber approximations for a slice of galaxies at spec_z of width dz.                                        
    eqn. (5) of https://arxiv.org/pdf/1511.04457.pdf                                                                                                      
    '''
    ks        = (Llls + 0.5)/comoving_distance(spec_z)                                              ## For the Phh evaluation in the integral, we take a zeff approx.
                                                                                                    ## i.e. \int dz .... Phh(zeff).                              
    result    = np.broadcast_to(Pab(ks, ba, bb, alpha, sigma, nowiggle), (len(z), len(Llls)))       ## Broadcast to each redshift.                                     
  else:
    result    = np.zeros((len(spec_z), len(Llls)))

    for i, j in enumerate(z):
      ks      = (Llls + 0.5)/chis[i]                                                                ## Evaluate Phh for each z and k(ell, z).               
                                                                                                    ## Given L mixes at range of k -- washes out BAO.  

      ## accounts for spatial to angular mapping as function of z, i.e. k = (l + 0.5)/chi(z) but neglects redshift evolution of Pps(k).
      result[i,:] = Pab(ks, ba, bb, alpha, sigma, nowiggle)                                         

  prefactor   = (cosmo.H(z).value/const.c.to('km/s').value)*(sliced_pz(z, spec_z, survey_dzs[surveya], surveya, True)/chis) \
                                                           *(sliced_pz(z, spec_z, survey_dzs[surveyb], surveyb, True)/chis)

  integrand   =  prefactor[:, None]*result
  integrand  /=  params['h_100']                                                                    ## account for [h^-1 Mpc]^3 of P(k) and h^-1 Mpc of chi_g.
 
  return  simps(integrand, dx = z[1] - z[0], axis=0)                                                ## integral over z. 

def _aCab(alpha, Llls, spec_z, surveya="LSST", surveyb="QSO", sigma=2., zeff=True, nowiggle=True):
    ## wrapper with alpha as leading argument for derivative calc.
    return Cab(Llls, spec_z, surveya, surveyb, alpha, sigma, zeff=True, nowiggle=nowiggle)

def _sCab(sigma, Llls, spec_z, surveya="LSST", surveyb="QSO", alpha=1., zeff=True, nowiggle=True):
    ## wrapper with sigma as leading argument for derivative calc.                                                                                         
    return Cab(Llls, spec_z, surveya, surveyb, alpha, sigma, zeff=True, nowiggle=nowiggle)

def get_errors(Fisher, nowiggle=False, print_it=True):
    iFisher   = inv(Fisher)
    nowiggle  = str(nowiggle)

    if  Fisher[0][0] > 0.0:
        CondErr  = 1./np.sqrt(Fisher[0][0])
    else:
        CondErr  = np.NaN 

    if iFisher[0][0] > 0.0:
        MargErr  =    np.sqrt(iFisher[0][0])
    else:
        MargErr  = np.NaN 

    if print_it == True:
        print "\n\nFisher matrix for no_wiggle = %s:  " % nowiggle
        print  Fisher

        print "\nInverse Fisher:"
        print  iFisher
        
        print "\nConditional error on alpha: %.6lf" % CondErr
        print "Marginal    error on alpha: %.6lf"   % MargErr

    return iFisher, CondErr, MargErr


if __name__ == "__main__":
    surveys  = "CMASS"
    surveyp  = "SDSS9"

    spec_z   = 0.75
    fsky     = 0.23                                                    ## Assumes CMASS is the limiting area.   
                                                                       ## QSOs are slightly smaller, but same ball park.
    alpha    = 1.0

    sigma_z0 = 7.0                                                     ## [Mpc/h] at z=0 for LCDM models (close to what we now believe).
    sigma    = sigma_z0 * growth_factor(1./(1. + spec_z))              ## Scales as D(z).

    ells     = np.arange(2., 1500., 1.)

    results  = {'False': {}, 'True': {}, 'Diff.': {}} 

    pl.clf()

    print "\nPerpendicular displacement at z = %.3lf:  %.4lf [Mpc/h]\n" % (spec_z, sigma)

    ## types = [False]
    types = [False, True]

    for nowiggle in types:
        Css      = Cab(ells, spec_z, surveys, surveys, alpha, sigma, zeff=True, nowiggle = nowiggle)
        Cps      = Cab(ells, spec_z, surveyp, surveys, alpha, sigma, zeff=True, nowiggle = nowiggle)
        Cpp      = Cab(ells, spec_z, surveyp, surveyp, alpha, sigma, zeff=True, nowiggle = nowiggle)

        Nss      = Ngg(ells, spec_z, survey_dzs[surveys], surveys)
        Npp      = Ngg(ells, spec_z, survey_dzs[surveyp], surveyp)           ## Read in dz from surveys dict. 

        var      = ((Cpp + Npp) * (Css + Nss) + Cps**2.) / (2. * ells + 1.)  
        ## var   = ( Cpp        *  Css        + Cps**2.) / (2. * ells + 1.)  ## Sample variance limited. 

        var     /= fsky

        err      = np.sqrt(var)

        if nowiggle == False:
            pl.semilogy(ells, prefactor(ells) * Css, 'k-',  label=r'$C_{ss}$')
            pl.semilogy(ells, prefactor(ells) * Cpp, 'r-',  label=r'$C_{pp}$') 
            pl.semilogy(ells, prefactor(ells) * Cps, 'b-',  label=r'$C_{ps}$')

            pl.semilogy(ells, prefactor(ells) * Nss, 'k--', label=r'$N_{ss}$')
            pl.semilogy(ells, prefactor(ells) * Npp, 'r--', label=r'$N_{pp}$')

            pl.semilogy(ells, prefactor(ells) * err, 'b--', label=r'$\sigma_{ps}$')
        
            pl.xlabel(r'$L$')
            pl.ylabel(r'$(2L+1) \ C_{\rm{ps}} (L) \ / \ 4 \pi$')
    
            pl.legend(loc=4, ncol=2)
            
            pl.title("%s x %s @ z = %.2lf" % (surveyp, surveys, spec_z))
        
            pl.savefig('Cps.pdf')

        ## ** Fisher matrix calculation. **
        ## n:   order of derivative; order: e.g. 5-point derivative, must be odd.                                                      
        ## see: https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.misc.derivative.html
        dCldA   = derivative(_aCab, alpha, dx=1e-6, n=1, order=5, args=(ells, spec_z, surveyp, surveys, sigma, True, nowiggle))
        dCldS   = derivative(_sCab, sigma, dx=1e-6, n=1, order=5, args=(ells, spec_z, surveyp, surveys, alpha, True, nowiggle))
    
        FAA     = np.sum(dCldA * dCldA / var)
        FAS     = np.sum(dCldA * dCldS / var)
        FSS     = np.sum(dCldS * dCldS / var)
        
        Fisher  = np.array([[FAA, FAS], [FAS, FSS]])

        ## Patej & Eisenstein prior:  0.8 < \alpha < 1.25
        prior_width = 0.2
        
        ## In this case, if \sigma is the width of the prior pdf, simply add 1/\sigma^2 to 
        ## the on-diagonal element corresponding to that variable (https://arxiv.org/pdf/0906.4123.pdf).
        Fisher[0][0]              += 1./prior_width**2.

        iFisher, CondErr, MargErr  = get_errors(Fisher, nowiggle)
     
        results[str(nowiggle)]     = {'Fisher': Fisher, 'iFisher': iFisher, 'CondErr': CondErr, 'MargErr': MargErr}

    results['Diff.']['Fisher']     = results['False']['Fisher']
    results['Diff.']['Fisher']    -= results['True']['Fisher']

    get_errors(results['Diff.']['Fisher'], "Diff.")

    print "\n\n***  Patej and Eisenstein  ***"
    print "D_M(z = 0.64) = (2418 \pm 73 Mpc) (rs / rs, fid)"
    print "\nD_M frac. err.: %.6lf" % (73./2418.)

    print "\nDone\n"
