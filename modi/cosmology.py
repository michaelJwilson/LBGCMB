## This class takes in a cosmology from astropy and H0 value in order to define some "extra" functions for that cosmology. 
## If H0 is none, take from the astropy-cosmology. Specify if you want to use 100. 
import  numpy
import  sys, os

from    scipy.integrate   import romberg, quad
from    scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from    astropy.cosmology import FlatLambdaCDM

from    modi_halofit      import halofit 
import  modi_tools        as     tools

pwd  = os.getcwd()

sys.path.append(pwd)
 
class Cosmology():
    ## To use scale factor, explicitly specify a = xx. Otherwise, redshift is assumed.
    def __init__(self, M = 0.292, L = None, H0 = 100., B = None, sig8 = 0.819, klin = None, plin = None, pfile = None):
        self.M = M

        if L is None: ## Flat LCDM universe.
            L = 1 - M

        self.L     = L
        self.B     = B  ## om_b
        self.H0    = H0
        self.sig8  = sig8
        self.cin   = 1000./(299792458.)  ## 1 / speed of light. 

        self.cosmo = FlatLambdaCDM(H0 = H0, Om0 = M, Ob0 = B, name = "cosmo")

        self.pfile = pfile

        if pfile is not None:
            self.klin, self.plin = numpy.loadtxt(pfile, unpack = True)

        else:
            self.klin = klin  ## Pk file not provided, expected on call. 
            self.plin = plin
        
        self.Dplus = self._eval_Dplus()
    
    def _za(self, z, a):
        if a is None:
            if z is None:
                print("At least one of scale factor or redshift is needed.")

                return None
            else: 
                a = self.ztoa(z)

                return z, a
        else:
            z = self.atoz(a)

            return z, a
    
    def atoz(self, a):
        return 1./a - 1.

    def ztoa(self, z):
        return 1./(1. + z)

    def OmegaMa(self, z = None, a = None):
        """ return Omega_m(a) """
        M, L = self.M, self.L
        z, a = self._za(z, a)

        return M/a**3./(M/a**3. + L)

    def Fomega1(self, a):
        """ Return \omega_m**(3./5.) """
        M  = self.M
        L  = self.L
        H0 = self.H0

        omega_a = M / (M + (1. - M - L) * a + L * a ** 3.)
        
        return omega_a**(3./ 5.)

    def Fomega2(self, a):
        """ Return 2. * \omega_m**(4./7.) """
        M  = self.M
        L  = self.L
        H0 = self.H0

        omega_a = M / (M + (1. - M - L) * a + L * a ** 3.)

        return 2. * omega_a ** (4./ 7.)

    def _eval_Dplus(self):
        """ Return un-normalized D(a).... evaluated once because lazy"""
        M       = self.M
        L       = self.L
        H0      = self.H0

        logamin = -20.
        Np      = 1000

        logx    = numpy.linspace(logamin, 0, Np)  ## Can use np.logspace with base. 
        x       = numpy.exp(logx)

        def kernel(loga):
            a = numpy.exp(loga)
            
            return (a * self.Ea(a = a)) ** -3. * a ## da = a * d loga

        y = self.Ea(a = x) * numpy.array([romberg(kernel, logx.min(), loga, vec_func=True) for loga in logx])

        return interpolate(x, y)

    def Dgrow(self, z = None, a = None):
        """ return normed: D(a)/D(a = 1.) """
        z, a = self._za(z, a)

        return self.Dplus(a)/self.Dplus(1.)

    def Ea(self, z = None, a = None):
        """ Return H(a)/H0 """
        z, a = self._za(z, a)

        M, L = self.M, self.L

        return a ** -1.5 * (M + (1. - M - L) * a + L * a ** 3.) ** 0.5

    def Ha(self, z = None, a = None):
        """ Return H(a) """
        z, a = self._za(z, a)

        return self.H0*self.Ea(a = a)

    ## Comoving distance, chi(a). 
    def chia(self, z = None, a = None):
        """ Return comoving distance to scale factor a in Mpc/h"""
        fac   = 100.*self.cin

        f     = lambda x: (self.Ea(a = x)*x**2)**-1.
        z, a  = self._za(z, a)

        if type(a) == numpy.ndarray:
            y       = numpy.zeros_like(a)

            for foo in range(len(a)):
                y[foo] = quad(f, a[foo], 1)[0]

            return y/fac

        else:
            return quad(f, a, 1)[0]/fac

    def k_chil(self, l, z = None, a = None):
        ''' k corresponding to \ell at distance chi(z) '''
        z, a = self._za(z, a)

        return (l + 0.5)/self.chia(a = a)  ## Added first order Limber correction, (ell + 0.5)/chi; no longer provided in call - argh!  

    def pkalin(self, z = None, a = None):
        ''' Calculate linear power spectrum at a(z) '''
        z, a = self._za(z, a)        

        return self.klin, self.plin*self.Dgrow(a = a)**2.

    def pkanlin(self, z = None, a = None, sig8 = None):
        ''' Calculate non-linear power spectrum at a, assuming halofit form '''
        z, a = self._za(z, a)        

        if sig8 is None:
            sig8 = self.sig8

        if z > 4:
            return self.pkalin(z, a)

        else:
            pdiml = self.plin*self.klin**3/(2.*numpy.pi**2)
            pgive = pdiml*self.Dgrow(a = a)**2.
        
            ## Halofit returns Del^2.
            pndim = halofit(k = self.klin, delta_k = pgive, sigma_8 = sig8, z = self.atoz(a), cosmo = self.cosmo, takahashi = True)
            pnlin = pndim * (2.*numpy.pi**2)/self.klin**3
        
            return self.klin, pnlin
                
    def ppota(self, z = None, a = None,  nlin = False):
        ''' Calculate the power spectrum of the potential at a '''
        z, a = self._za(z, a)
        
        if nlin: 
            pmat =  self.pkanlin(a = a)[1]

        else:
            pmat =  self.pkalin(a = a)[1]

        return self.klin, (9 * self.M**2 * self.H0**4 * pmat)/(4 * a**2 * self.klin**4) * self.cin**4
