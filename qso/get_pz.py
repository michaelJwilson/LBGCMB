import  os
import  pickle
import  numpy   as      np
import  pylab   as      pl

from    utils   import  fwhm


root  = os.environ['LBGCMB']
 
def load_survey(survey="CMASS"):
    import  astropy.io.fits  as  fits

    """
    Load all redshifts for a given catalogues. 
    Currently available:  CMASS, QSO and SDSS9.
    """
    all_zs = []

    print("\n\nLoading survey %s" % survey)

    if survey == "CMASS":
        for area in ["North", "South"]:
            data  = fits.open(root + '/cmass/galaxy_DR12v5_CMASS_%s.fits' % area)
            
            zs    = data[1].data['Z']
            wsys  = data[1].data['weight_systot']         ## Star and seeing corrections.  
            wnoz  = data[1].data['weight_noz']            ## https://data.sdss.org/datamodel/files/BOSS_LSS_REDUX/galaxy_DRX_SAMPLE_NS.html
            wcp   = data[1].data['weight_cp']             ## https://arxiv.org/pdf/1203.6594.pdf
            
            wght  = wsys * (wnoz + wcp - 1.0)

            for zee in zs:
                all_zs.append(zee)
            
    elif survey == "QSO":
        data    = fits.open(root + '/qso/DR12Q.fits')

        for z in data[1].data:
            all_zs.append(z[8])
    
    elif survey == "SDSS9":
        for fname in ['R+00_+01.fits', 'R+01_+02.fits', 'R+34_+35.fits', 'R+35_+36.fits']:
            data       = fits.open(root + '/sdss9/' + fname)
            zs         = data[1].data['zphot']
        
            for zee in zs:
              all_zs.append(zee)
    
    else:
        raise  ValueError("\n\nProvided survey is not available.")

    all_zs = np.array(all_zs)

    calc_pz(all_zs, survey=survey)

    return all_zs

def calc_pz(zs, survey=None):
    """
    Given an array of redshifts, calculate p(z).
    If survey is defined, write to a pickle file.
    """

    bins          =   np.arange(0.0, 6.0, 0.01)
    (dNdz, bins)  = np.histogram(zs, bins=bins)

    dNdz          = dNdz.astype(np.float)
    bins          = bins[:-1]

    dz            = bins[1] - bins[0]
    midz          = bins + dz/2.

    dNdz         /= dNdz.sum()
    pz            = dNdz/dz                               ## sampled at each midz; assumed constant across the bin.                                        

    if survey is not None:
      print("\n\nPickling p(z) for %s" % survey)

      pickle.dump(np.column_stack((midz, pz)), open(survey + '/dNdz.p', 'wb'))

def get_pz(z, survey = "CMASS"):
    from  scipy.interpolate  import interp1d
    
    """
    Get p(z) for given survey; loading the relevant pickle file. 
    """

    surveys      = {"DESI-QSO": {"fname":   'qso/dNdz.p'},   "QSO": {"fname":   'qso/dNdz.p'},\
                       "SDSS9": {"fname": 'sdss9/dNdz.p'}, "CMASS": {"fname": 'CMASS/dNdz.p'}}

    fpath        =  root + surveys[survey]["fname"]

    data         =  pickle.load(open(fpath, 'r'))

    midz         =  data[:,0]
    pz           =  data[:,1]

    pz_interp    =  interp1d(midz, pz, bounds_error=False, fill_value=0.0)

    return  pz_interp(z)

def nbar(survey = "CMASS", printit=False):
  """
  Get the number of galaxies per sq. deg for given survey.
  """
  if survey == "CMASS":
    N            = 618806. + 230831.                  ## Observed CMASS numbers -- sum of north and south galactic caps.                            
    area         =             9493.                  ## deg^2                                                                                     
  
  elif survey == "QSO":
    N            = 297301.                            ## Observed quasar numbers;                                                                 
    area         = 9376.                              ## deg^2                                                                                         

  elif survey == "SDSS9":
    N            = 6594677.                           ## Observed SDSS9 galaxies (for Patej and Eisenstein.);                                          
    area         = 14555.                             ## deg^2                                                                                        

  elif survey == "DESI-QSO":
    N            = 2.4e6
    area         = 14000.

  else:
    raise  ValueError("Spec. x Photo. p(z) normalisation is not defined for %s." % survey)

  if printit:
      print "\n%s survey has %.6lf g/deg2." % (survey, N / area)

  return  N / area                                    ## Number of galaxies per sq. degree.
 
def zsplit_QSO_samplestats(zs):
    area      = 9376.                                 ## deg^2  
    
    NTo       = len(zs)
    
    Nlo       = len(zs[zs < 2.0])
    Nhi       = len(zs[zs > 2.0])
    
    lonbar    = Nlo / area
    hinbar    = Nhi / area

    mean_zlo  = np.sum(zs[zs < 2.0]) / Nlo
    mean_zhi  = np.sum(zs[zs > 2.0]) / Nhi

    std_zlo   = np.std(zs[zs < 2.0])
    std_zhi   = np.std(zs[zs > 2.0])

    fwhm_zlo  = fwhm(std_zlo)
    fwhm_zhi  = fwhm(std_zhi)

    print lonbar, hinbar
    print mean_zlo, mean_zhi
    print std_zlo, std_zhi
    print fwhm_zlo, fwhm_zhi


if __name__ == "__main__":
    print "\nWelcome to a (Patej & Eisenstein) spec. x photo z calculator.\n\n"
    
    dz   = 0.01    
    zs   = np.arange(0.0, 6.0, dz)

    ## Also available: SDSS9, DESI-QSO.
    for survey, zcut in zip(["QSO"], [3.0]):
    ## for survey, zcut in zip(["QSO", "CMASS"], [3.0, 0.6, None, None]):
        zs = load_survey(survey)

        zsplit_QSO_samplestats(zs)

        '''
        nbar(survey, printit=True)

        pz = get_pz(zs, survey)

        for i, zee in enumerate(zs):
            print zee, pz[i]

        if zcut is not None:
            print "\n\np(z) estimate of %s with z > %.1lf:  %.4lf\n" % (survey, zcut, dz * pz[zs > zcut].sum())
        '''
    print("\n\nDone.\n\n")
