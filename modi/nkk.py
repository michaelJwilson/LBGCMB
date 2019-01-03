# Script to calculate the instrumental and Fisher Noise for TT & EB
# for a given CLASS data, and \Delta_T and \Beam-Size
# You can specify --temp & --beam on command line
# For input files and output folder, easiest is to edit in class Files
# However you can also pass these through command line

import  numpy
import  pylab              as pl
import  sys, argparse

from    scipy.integrate    import simps
from    scipy.interpolate  import InterpolatedUnivariateSpline as interpolate
from    DetectorNoise      import DetectorNoise


def modi_lines(fname = "runpb2_lsst_z200_bw_15_dT_10.txt"):
  cols = ['l', 'Clkk', 'Nkk', 'Clgg', 'Clgg_zeff', 'shot', 'Clkg', 'Clkg_zeff', 'var_kg', 'var_kg_zeff', '']
  data = pd.read_csv(fname, header=None, names=cols, delim_whitespace=True)

  pl.plot(data['l'].values, prefactor(data['l'].values)*data['Nkk'], 'k--', label='Chirag')

def N_inst(thetab, delT, l, mode = 'T'):
    ''' Return N_l spectrum given Delta_T [uK-arcmin] and beam-size [arcmin]'''
    arcm    =  (60*180/numpy.pi)**-1
    thetab *= arcm

    Tcmb    = 2.786e6
    Nl      = (delT*arcm/Tcmb)**2 * numpy.exp( l*(l+ 1)*thetab**2. /(8* numpy.log(2)))

    if mode == 'T':
        return Nl
    else:
        return 2.*Nl

def N_fisher(L, l, icl1, icl2, icln1, icln2, mode):
    '''Do the variance double integral for TT, EE and EB'''
    # Will break Simpson's rule.
    #l      = l[l!=L] #  remove the same l=L since then l2 is undefined; i.e. L = 0, for which the lesing potential, \psi(0), is unobservable - no deflection.
    nphi   = 200
    
    phis   = numpy.linspace(0, 2*numpy.pi, nphi)
    cosphi = numpy.cos(phis)

    cl1    =  icl1(l)  # interp Cl(ell)
    cln1   = icln1(l)  # noise added

    # Generate l2 \equiv |l - L|
    Llcosphi =  L*l.reshape(-1, 1)*cosphi.reshape(1, -1)         # Col vector:      L*l * cos(phi)                -> outer product
    l2       = (L**2 + (l**2).reshape(-1, 1) - 2*Llcosphi)**0.5  # Col vector: sqrt(L**2. + l**2. - 2*l*cos(phi)) -> outer product

    try:
        cl2 = icl2(l2)

    except:
        cl2 = 0.0                                               # Zero outside interpolation range.
    
    cln2    = icln2(l2)                                         # No exception for noise; added power at high l2. 

    # Evaluate the fraction at the l1, l2 grid. Mode dependent
    num = (Llcosphi*cl1.reshape(-1, 1) + (L**2 - Llcosphi)*cl2)**2.  # col vector: {[L**2. - L*l*cos(phi)]*Cl(l2) + L*l*cos(phi)*Cl(l1)}**2. 
    den = cln1.reshape(-1, 1)*cln2                                   # l1 -> l - L, l2 = l    

    ## TEST
    #num = (Llcosphi*cl1.reshape(-1, 1))**2.                         # col vector: {[L**2. - L*l*cos(phi)]*Cl(l2) + L*l*cos(phi)*Cl(l1)}**2.           
    #den = cln1.reshape(-1, 1)*cln2

    if mode == "TT":
      den   *= 2.
    
    else:
        # Angles between 2 l's
        sinphi12 =  L*numpy.sin(phis).reshape(1, -1)/l2            # sin(phi_12) =  L*sin(phi)/l2 
        cosphi12 = (L*cosphi.reshape(1, -1) - l.reshape(-1, 1))/l2 # cos(phi_12) = [L*cos(phi) - l]/l2

        if mode == "EB":
            num *= (2*sinphi12*cosphi12)**2.

        elif mode == "EE":
            num *= (cosphi12**2. - sinphi12**2.)**2.
            den *= 2.

    frac = num/den         

    frac[l2 > l.max()] = 0.
    frac[l2 < l.min()] = 0.

    # Integrate along phi and then ell.
    int1 = l * simps(frac, phis, axis = -1)

    int2 = simps(int1, l)

    return int2

def N_fisherTE(L, l, iclte, iclnt, iclne, iclnte):
    '''Do the variance double integral for TE'''
    # Will break Simpson's rule.
    #l        = l[l!=L]                              # Remove the same l=L since then l2 is undefined

    nphi      = 200

    phis      = numpy.linspace(0, 2*numpy.pi, nphi) 
    cosphi    = numpy.cos(phis)

    Llcosphi  =  L*l.reshape(-1, 1)*cosphi.reshape(1, -1)            #  L * lx 
    l2        = (L**2. + (l**2.).reshape(-1, 1) - 2.*Llcosphi)**0.5  # |l - L|

    Ll2cosphi =  L**2. - Llcosphi                                    # L^2 - L * lx

    sinphi12  =  L*numpy.sin(phis).reshape(1, -1)/l2
    cosphi12  = (L*cosphi.reshape(1, -1) - l.reshape(-1, 1))/l2

    cos2phi12 = cosphi12**2. - sinphi12**2.

    #  For F (numerator)
    cl1te  = iclte(l).reshape(-1, 1)  # unlensed Cl_TE(l)
    cl2te  = iclte(l2)                # unlensed Cl_TE(|l-L|) 

    cl1ntt = iclnt(l).reshape(-1, 1)  # lensed Cl_TT with noise at l
    cl2ntt = iclnt(l2)                # lensed Cl_TT with noise at |l-L|

    cl1nee = iclne(l).reshape(-1, 1)  # lensed Cl_EE with noise at l 
    cl2nee = iclne(l2)                # lensed Cl_EE with noise at |l-L|

    cl1nte = iclnte(l).reshape(-1, 1) # lensed Cl_TE with NO noise at l
    cl2nte = iclnte(l2)               # lensed Cl_TE with NO noise at |l-L|
                                      
    fl1l2  = Llcosphi * cl1te * cos2phi12 + Ll2cosphi*cl2te  # L * lx * Cl_TE(l) * cos(2 * phi_12) + (L^2 - L * lx) * Cl_TE(l - L)
    fl2l1  = Llcosphi * cl1te + Ll2cosphi * cl2te * cos2phi12  # L * lx * Cl_TE(l) * cos(2 * phi_12) + (L^2 - L * lx) * Cl_TE(l - L)

    num    = cl1nee * cl2ntt * fl1l2 - cl1nte * cl2nte *fl2l1             # Cl_EE(l) with noise * Cl_TT(|l-L|) with noise * fl1l2 - Cl_TE(l) with noise
    den    = cl1ntt * cl2nee * cl1nee * cl2ntt - (cl1nte * cl2nte)**2.   

    F      = num/den

    frac   = fl1l2*F

    # Set to zero outside range
    frac[l2 > l.max()] = 0.
    frac[l2 < l.min()] = 0.

    # Integrate along the angles and then "l"
    int1 = l * simps(frac, phis, axis = -1)
    
    int2 = simps(int1, l)   
     
    return int2

class CMB:
    def __init__(self, clfile="RunPB02_cl.dat", clsfile="RunPB02_cl_lensed.dat"):
        self.clfile  = clfile
        self.clsfile = clsfile
        
        self.read_class()
        ## self.set_camb()
        
        # Value of ells at which to calculate noise, small array for debug run
        self.lint   = numpy.linspace(1., 3000., 10000)
        self.lintP  = numpy.linspace(1., 5000., 10000)
        
        self.lat    = numpy.logspace(2.0,  3.5,    15) # numpy.logspace(1, 3.47,   500)

        self.kappa  =  self.cl[:, 5]*self.ll**2 /4.
        self.kappas = self.cls[:, 5]*self.lls**2 /4.

    ## def set_camb(self):
    ##  (self.camb_lCls, self.camb_uCls) = prep_Clxy() 

    def read_class(self):
        '''Read class output file for unlensed(clfile) & lensed(clsfile) data.
        The assumed format of columns is:  l, TT, EE, TE, BB, phiphi, TPhi, Ephi
        '''
        # Columns
        self.cols = ['l','TT', 'EE', 'TE', 'BB', 'phiphi', 'TPhi', 'Ephi' ]
        self.coldict = dict()

        for foo in range(len(self.cols)):
            self.coldict[self.cols[foo]] = foo

        # Read files
        self.cl  = numpy.loadtxt(self.clfile) 
        self.cls = numpy.loadtxt(self.clsfile) 

        # Convenient variables
        self.ell  = self.cl[:, 0]
        self.ells = self.cls[:, 0] #  Class outputs unlensed data for 500 more 'l'.
        self.ll   = self.ell*(self.ell   + 1)
        self.lls  = self.ells*(self.ells + 1)
        self.tpi  = 2*numpy.pi

        # Remove factors of ell from class-output
        self.cl[:, 1:]  =  self.cl[:, 1:] *(self.tpi/self.ll.reshape(-1, 1))
        self.cls[:, 1:] = self.cls[:, 1:] *(self.tpi/self.lls.reshape(-1, 1))

    def do_noise(self, mode1, mode2, thetab, delT):
        ''' Return Fisher noise for mode = mode1 + mode2 '''
        #  Signal-spectra in unlensed autospectrum for component modes
        auto1 = mode1 + mode1
        auto2 = mode2 + mode2

        if mode1 + mode2 == "TE" or mode1 + mode2 == "ET":
            print("Wrong function")

            return None

        icl1   = interpolate(self.ell, self.cl[:, self.coldict[auto1]])
        icl2   = interpolate(self.ell, self.cl[:, self.coldict[auto2]])
        
        # Noise-spectra in lensed autospectrum for component modes + Instrumental Noise (2*Nl for E&B)
        Nl1    = DetectorNoise(self.ells, thetab, delT, auto1)  # N_inst(thetab, delT, self.ells, mode1)
        Nl2    = DetectorNoise(self.ells, thetab, delT, auto2)  # N_inst(thetab, delT, self.ells, mode2)

        icln1  = interpolate(self.ells, self.cls[:, self.coldict[auto1]] + Nl1)
        icln2  = interpolate(self.ells, self.cls[:, self.coldict[auto2]] + Nl2)

        noise  = numpy.zeros_like(self.lat)

        if mode1 + mode2 == "EB"  or mode1 + mode2 == "EE":
            print('Go to higher lmax for %s' % (mode1 + mode2))

            for foo in range(len(noise)):           
                noise[foo] = N_fisher(self.lat[foo], self.lintP, icl1, icl2, icln1, icln2, mode1 + mode2)
        else:
            for foo in range(len(noise)):
                noise[foo] = N_fisher(self.lat[foo], self.lint, icl1, icl2, icln1, icln2, mode1 + mode2)

        noise = (2*numpy.pi)**2./noise 

        return Nl1, Nl2, noise * (self.lat*(self.lat + 1))**2. /4.


    #  Keep TE separate
    def do_noiseTE(self,  thetab, delT):
        '''Return Fisher noise for mode = 'mode1mode2'
        '''
        #  Signal-spectra in unlensed autospectrum for component modes
        mode1  = "T"
        mode2  = "E"
        
        icl3   = interpolate(self.ell, self.cl[:, self.coldict["TE"]])
        #icl3   = self.camb_uCls['TE']
        
        #  Noise-spectra in lensed autospectrum for component modes + Instrumental Noise (2*Nl for E&B)
        Nl1    = DetectorNoise(self.ells, thetab, delT, mode1 + mode1) # N_inst(thetab, delT, self.ells, mode1)
        Nl2    = DetectorNoise(self.ells, thetab, delT, mode2 + mode2) # N_inst(thetab, delT, self.ells, mode2)
        
        icln1  = interpolate(self.ells, self.cls[:, self.coldict["TT"]] + Nl1)
        icln2  = interpolate(self.ells, self.cls[:, self.coldict["EE"]] + Nl2)
        icln3  = interpolate(self.ells, self.cls[:, self.coldict["TE"]])          # NO noise added to TE.

        #icln1  = self.camb_lCls['TT']
        #icln2  = self.camb_lCls['EE']
        #icln3  = self.camb_lCls['TE']        # NO noise added to TE.

        noise  = numpy.zeros_like(self.lat) # noise for each L (log spacing).  
        
        for foo in range(len(noise)):
            noise[foo] = N_fisherTE(self.lat[foo], self.lint, icl3, icln1, icln2, icln3)

        noise = (2*numpy.pi)**2./noise

        return Nl1, Nl2, noise * (self.lat*(self.lat + 1))**2. / 4.
       
    def plot_Cls(self, mode1, mode2, thetab, delT):
          Nl1    = DetectorNoise(self.ells, thetab, delT, auto1) # N_inst(thetab, delT, self.ells, mode1)
          Nl2    = DetectorNoise(self.ells, thetab, delT, auto2) # N_inst(thetab, delT, self.ells, mode2)

          icl1   = interpolate(self.ell, self.cl[:, self.coldict[auto1]])
          icl2   = interpolate(self.ell, self.cl[:, self.coldict[auto2]])

          icln1  = interpolate(self.ells, self.cls[:, self.coldict[auto1]] + Nl1)
          icln2  = interpolate(self.ells, self.cls[:, self.coldict[auto2]] + Nl2)

          ells   = numpy.logspace(0.0, 4., 5000)

          pl.loglog(ells, prefactor(ells, 2)*icl1(ells),  'm--', label='Modi')
          pl.loglog(ells, prefactor(ells, 2)*icln1(ells), 'm--', label='Modi')

    def plot_TECls(self, mode1, mode2, thetab, delT):
          Nl1    = DetectorNoise(self.ells, thetab, delT, 'TE') # N_inst(thetab, delT, self.ells, mode1)                                               

          icl1   = interpolate(self.ell, self.cl[:,   self.coldict['TE']])
          icln1  = interpolate(self.ells, self.cls[:, self.coldict['TE']] + Nl1)

          ells   = numpy.logspace(0.0, 4., 5000)

          pl.plot(ells, prefactor(ells, 2)*icl1(ells),  'm--', label='Modi')
          pl.plot(ells, prefactor(ells, 2)*icln1(ells), 'c--', label='Modi')


def modi_calc(x, xp):
  # Callable: TT, EE, EB, TE; only TE requires cross spectra. 

  cmb = CMB()

  if x+xp == 'TE':
    Ninst_x, Ninst_xp, Nfish_xxp = cmb.do_noiseTE(1., 1.)
  
  else:
    Ninst_x, Ninst_xp, Nfish_xxp = cmb.do_noise(x, xp, 1., 1.)
  
  pl.semilogy(cmb.lat, prefactor(cmb.lat)*Nfish_xxp, 'k--')
