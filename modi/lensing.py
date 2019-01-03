import  math
import  numpy              as np
import  pylab              as pl
import  sys, os, argparse

import  modi_pk_pregen     as PK
import  modi_cosmology     as cosmology 
import  modi_halofit       as halofit
import  modi_tools         as tools 

from    scipy.integrate    import quad, simps
from    scipy.interpolate  import InterpolatedUnivariateSpline as interpolate


pwd = os.getcwd()
sys.path.append(pwd)

def prefactor(ells, n=1):
  if n == 0:
    return 1.0

  if n == 1:
    return (2.*ells + 1.0)/(4.*np.pi)

  if n == 2:
    return ells*(ells + 1.0)/(2.*np.pi)

## If file paths are not specified in the command line arguments, use these:
class Files:
    parentdir  = "/Users/M.J.Wilson/work/cmbcross"
    pkfile     = parentdir + "/data/RunPB/cambPB_pk.txt"  ## Linear P(k).
    outpath    = parentdir + "/data/RunPB/"

    classfile  = None
    errorfile  = None

class CosmoPar:
    M  = 0.292
    L  = 1 - M
    H0 = 100.
    h  = 0.69

class LensingX:
    ''' Lensing cross-spectra in the thin-slice approximation.'''
    def __init__(self, z, dndz, cosmo = None, l = None, dz = 0.5, nz = 100, ggnorm = None):
        if cosmo is None:
            cosmo = cosmology.Cosmology(M = CosmoPar().M, pfile = Files().pkfile)
        
        if l is None:
            l     = np.logspace(1, 3.47, 500) ## np.linspace(100, 4000, 200)
        
        if dndz is None:
            dndz  = lambda x: 1
        
        self.ggnorm = ggnorm                             ## Integral constraint; stepwise approx to dN/dz. 
    
        self.z      = z
        self.l      = l 
        self.dz     = dz
        self.nz     = nz
        self.kknz   = 20*self.nz
        self.dndz   = dndz
        self.cosmo  = cosmo

        self._clsetup()

    def _clsetup(self):        
        zmin        = self.z - self.dz/2.
        zmax        = self.z + self.dz/2.

        self.zs     = np.linspace(zmin,   zmax, self.nz)    ## zs spanning z - dz/2 to z + dz/2 in nz bins. 
        self.kkzs   = np.linspace(0.001, 1100., self.kknz) 

        self.chis   = self.cosmo.chia(z = self.zs)          ## comoving distance at z.    
        self.kkchis = self.cosmo.chia(z = self.kkzs)

        ## Kernels
        self.kerg   = self.kernel_g(z = self.zs, dndz = self.dndz)
        self.norm   = np.trapz(self.kerg, self.chis)
        
        if self.ggnorm is not None:
            ## Useful for photometric when vary dndz but keep total N constant
            print "External gg-norm given = %0.3f, gg-norm estimated from dndz = %0.3f; Overriding estimated norm. Be wary!!!" % (self.ggnorm, self.norm)

            self.norm = self.ggnorm
        
        self.kerg    /= self.norm                               ## Normalise W^g(\chi) either externally or internally. 

        self.kercmb   = self.kernel_cmb(z = self.zs)            ## W_cmb; eqn. 2.3 of Modi et al.  
                                                               
        self.ggwts    = self.kerg**2./self.chis**2.             ## W_g   * W_g /chi^2; eqn. 2.3 of Modi et al. 
        self.kgwts    = self.kercmb*self.kerg/self.chis**2.     ## W_cmb * W_g /chi^2; eqn. 2.3 of Modi et al.        
        self.kkwts    = self.kernel_cmb(z = self.kkzs)*self.kernel_cmb(z = self.kkzs)/self.cosmo.chia(z = self.kkzs)**2. ## Extended to last scatter.
        
        self.dggwtsdN = 2.*self.kerg/self.chis**2. * self.cosmo.Ha(z = self.zs)/self.norm   ## Fisher derivative wrt N(z_i).
        self.dkgwtsdN =    self.kercmb/self.chis**2. * self.cosmo.Ha(z = self.zs)/self.norm

        ## Integration mesh
        self.kmesh     = np.zeros([self.l.size,   self.zs.size])
        self.kk_kmesh  = np.zeros([self.l.size, self.kkzs.size])

        for foo in range(self.zs.size):
           ## For each z, calculate ks assuming Limber correction: k = (ell + 0.5)/chi(z); NOTE: Corrected from original implementation.  
           self.kmesh[:, foo]    = self.cosmo.k_chil(self.l, z = self.zs[foo])        

        for foo in range(self.kkzs.size):    
           self.kk_kmesh[:, foo] = self.cosmo.k_chil(self.l, z = self.kkzs[foo])

    def kernel_cmb(self, z = None, a = None):
        ## Returns W^\kappa, eqn. (2.4) of Modi et al
        z, a    = self.cosmo._za(z, a)      ## Given z or a, return the pair.
        chistar = self.cosmo.chia(z = 1100.) ## Redshift of cmb
        chi     = self.cosmo.chia(z)

        f       = chi*(chistar - chi)/(chistar)

        fac     = 3.*self.cosmo.M*self.cosmo.H0**2./(2.*a) 
        
        return f*fac*self.cosmo.cin**2.     ## cin is 1/c due to H0 (km/s/Mpc).

    def kernel_g(self, z, dndz):
        ## Returns W^g, eqn. (2.4) of Modi et al
        Ha = self.cosmo.Ha(z = z)           ## Returns H(a)
        
        return Ha*dndz(z)

    def cl_zeff(self, kp, auto = False):
        ''' Assumes the power spectrum is constant across the redshift bin '''
        pval   = np.interp(self.kmesh, *kp)       ## (k, P(k)) for e.g. hh, or hm; interpolate to k mesh (specified by ell and z).
                                                  ## Units of P(k), and distance?
        if auto:
            cl = pval*self.ggwts

        else:
            cl = pval*self.kgwts

        return np.trapz(cl, self.chis, axis = -1)  ## Integral over comoving distance. 

    def cl_zeff_bitracer(self, kp, dndz2):
        ''' Cross-spectra with two galaxy tracers; assumes constant power spectrum across redshift bin '''
        pval   = np.interp(self.kmesh, *kp)
        
        ## Kernels
        kerg2   = self.kernel_g(z =  self.zs, dndz = dndz2)
        kerg2  /= np.trapz(kerg2, self.chis)       ## No option for external normalisation.

        ggwts2  = self.kerg*kerg2/self.chis**2.
        cl      = pval * ggwts2

        return np.trapz(cl, self.chis, axis = -1)

    def clz(self, ikpz, auto = False):
        ''' Assuming interpolating function for Pk(z) across the redshift bin '''
        pval = np.zeros_like(self.kmesh)
        
        for foo in range(self.nz):
            kp           = ikpz(z = self.zs[foo])

            pval[:, foo] = np.interp(self.kmesh[:, foo], *kp)

        if auto:
            cl = pval * self.ggwts

        else:
            cl = pval * self.kgwts

        return np.trapz(cl, self.chis, axis = -1)

    def cl_kkz(self, ikpz):
        ''' Assuming interpolating function for Pk(z) across the redshift bin '''
        pval = np.zeros_like(self.kk_kmesh)
        
        for foo in range(self.kknz):
            kp           = ikpz(z = self.kkzs[foo])
            pval[:, foo] = np.interp(self.kk_kmesh[:, foo], *kp, left=0.0, right=0.0)

        cl = pval*self.kkwts

        return np.trapz(cl, self.cosmo.chia(z = self.kkzs), axis = -1)

def calc_runpb(z, dndz, lo):
    ''' Calculate spectra in the Run PB cosmology. '''
    iz        = z*100.
        
    parentdir = "/Users/M.J.Wilson/work/cmbcross"
    
    if lo:    ## low mass halo bin; 
        hh    = np.loadtxt(parentdir + '/data/RunPB/hh_fine_z%03d_lo.txt' % iz, unpack = True)
        hm    = np.loadtxt(parentdir + '/data/RunPB/hm_fine_z%03d_lo.txt' % iz, unpack = True)

    else:
        hh    = np.loadtxt(parentdir + '/data/RunPB/hh_fine_z%03d.txt'    % iz, unpack = True)
        hm    = np.loadtxt(parentdir + '/data/RunPB/hm_fine_z%03d.txt'    % iz, unpack = True)

    lenz      = LensingX(z = z, dndz = dndz)
    l         = lenz.l
    
    clggzeff  = lenz.cl_zeff(hh, auto =  True)
    clkgzeff  = lenz.cl_zeff(hm, auto = False)

    toret     = [l, clggzeff, clkgzeff]
    
    return toret
    
def calc_halofit(z, dndz):
    ''' Calculate spectra assuming Halofit. '''
    parentdir  = "/Users/M.J.Wilson/work/cmbcross"

    klin, plin = np.loadtxt(parentdir + "/data/RunPB/cambPB_pk.txt", unpack = True)
    
    cosmo      = cosmology.Cosmology(M = CosmoPar.M, L = CosmoPar.L, pfile = Files.pkfile) ## args.cosmo
    bz         = tools.BiasZ(cosmo).const(b = 1.438)

    ## mmz     = lambda z: [klin,  cosmo.pkalin(z = z)[1]]  ## Linear
    mmz        = lambda z: [klin, cosmo.pkanlin(z = z)[1]]  ## Non-linear

    hmz        = lambda z: [klin, cosmo.pkanlin(z = z)[1]*bz(z)]
    hhz        = lambda z: [klin, cosmo.pkanlin(z = z)[1]*bz(z)**2.]
    
    mm         = mmz(z)
    hh         = hhz(z)
    hm         = hmz(z)

    lenz       = LensingX(z = z, dndz = dndz)
            
    l          = lenz.l
            
    clkkz      = lenz.cl_kkz(mmz)
    
    clggz      = lenz.clz(hhz, auto =  True)  ## NOTE: Changed from hh.
    clkgz      = lenz.clz(hmz, auto = False)  ## NOTE: Changed from hm.
    
    clggzeff   = lenz.cl_zeff(hh, auto =  True)
    clkgzeff   = lenz.cl_zeff(hm, auto = False)
            
    toret      = [l, clkkz, clggz, clggzeff, clkgz, clkgzeff]
            
    return toret

def calc_lpt(z, dndz):
    ''' Calculate spectra assuming LPT. '''
    lpt   = PK.Pk8CLEFT(z = z)
    
    fits  = list(np.loadtxt('joint_fit_results.log')[2*z - 2, 1: -1])

    xfits = fits[:5]
    afits = fits[:3] + fits[-2:]

    hh    = lpt(afits, auto = True)
    hm    = lpt(xfits, auto = False)
    
    lenz  = LensingX(z = z, dndz = dndz)
    l     = lenz.l

    clggzeff = lenz.cl_zeff(hh, auto = True)
    clkgzeff = lenz.cl_zeff(hm, auto = False)

    toret    = [l,  clggzeff,  clkgzeff]
    
    return toret

def modi_lensingcalc():
    lo              = False

    ## args         = tools.parse_arguments(Files = Files, CosmoPar = CosmoPar)
    cosmo           = cosmology.Cosmology(M = CosmoPar.M, L = CosmoPar.L, pfile = Files.pkfile) ## args.cosmo
    dndz            = tools.DnDz().func("lsst")  ## tools.DnDz().func(args.survey)
    zg              = 1                          ## args.zg
    fsky            = 0.5
    beam            = 15.
    temp            = 10.
    survey          = "lsst"

    parentdir       = "/Users/M.J.Wilson/work/cmbcross/"

    ## ell, C_l^kappa, C_l^kappa(lensed file), N_fisher^TT, N_fisher^EE,  N_fisher^TE,  N_fisher^EB, N_Instr^T, N_Instr^E, N_Instr^B
    noisefile       = parentdir + "/data/RunPB/noise/noise_ext4_bw_%d_dT_%d.txt" % (beam, temp)
    noise           = np.loadtxt(noisefile).T

    ell             = noise[0]
    clkk            = noise[1]

    ## Nkk
    nkk             = (1./noise[3] + 1./noise[6])**-1. 
    
    if lo:
        suffix      = 'runpb-lo'
        b1          = np.loadtxt('/Users/M.J.Wilson/work/cmbcross/data/RunPB/joint_fit_results_lo.log')[2*zg - 2, 1] + 1

    else:
        suffix      = 'runpb2'
        b1          = np.loadtxt('/Users/M.J.Wilson/work/cmbcross/data/RunPB/joint_fit_results.log')[2*zg-2, 1] + 1

    print '\nb1 for %s = %.3lf\n' % (suffix, b1)

    ## Cls (RunPB)                                                                                                                                           
    ## l, clggzeff, clkgzeff = calc_runpb(z = zg, dndz = dndz, lo = lo)                                                                                    
    ## clggz, clkgz          = l*0, l*0  

    ## N_gg errors                                                                                                                                           
    nobj            = tools.nbar(zg - 0.25, zg + 0.25, survey = survey)
    arcm            = np.pi/60./180.
    shotnoise       = 1./(nobj/arcm**2.)

    ## Halofit
    l, clkkz, clggz, clggzeff, clkgz, clkgzeff = calc_halofit(z = zg, dndz = dndz)

    ## LPT
    ## l, clggzeff, clkgzeff = calc_lpt(z = zg, dndz = dndz)

    err_kgz         = np.sqrt(((clkkz + nkk)*(clggz + shotnoise) + (clkgz)**2.)/(2*l + 1)/fsky)

    ## pl.plot(l, prefactor(l)*clkkz,   'k--', label=r'$C_{kk}$')
    ## pl.plot(l, prefactor(l)*nkk,     'k--', label=r'')

    ## pl.plot(l, prefactor(l)*clggz,   'k--', label=r'$C_{gg}$')
    
    pl.plot(l, prefactor(l)*clkgz,   'k--', label=r'$C_{kg}$')
    pl.plot(l, prefactor(l)*err_kgz, 'k--', label=r'')

    pl.legend()

    pl.xscale('linear')
    pl.yscale('log')

    ## pl.savefig('modi_lensing.pdf')    
    

if __name__=="__main__":
    modi_lensingcalc()

