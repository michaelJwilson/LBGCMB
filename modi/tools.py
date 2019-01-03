import math
import argparse
import sys, os, argparse

import numpy             as np
import modi_cosmology    as cosmology
import modi_tools        as tools
import modi_pk_pregen    as PK

from   collections       import OrderedDict as odict

pwd  = os.getcwd()

sys.path.append(pwd)

def binit(l, y, lmin = 10, lmax = 2.1e3, dl = 0.1):
    lt    = lmin
    lbins = []

    while lt < lmax:
        lbins.append(lt)
        lt += lt*dl

    lbins = np.array(lbins)
    yt    = np.interp(lbins, l, y)

    return lbins, yt

def binl(lmin = 10, lmax = 2.1e3, dl = 0.1):
    lt    = lmin
    lbins = []
    
    while lt < lmax:
        lbins.append(lt)
        lt += lt*dl
    
    lbins   = np.array(lbins)
    
    return lbins

def sig8( k, p):
    kr = k*8
    wt = 3 * (np.sin(kr)/kr - np.cos(kr))/kr**2
    
    if wt is 0:
        wt = 1
    
    sig =  np.trapz(p * wt**2. * k**2, k)/2/np.pi**2
    
    return sig**0.5
    
def lagrange(z, zz):  
    c        = np.ones_like(zz)
    zzdist   = zz[:, None] - zz[None, :]

    singular = zzdist ==  0

    zzdist[singular]  = 1.0

    fac           = (z - zz)/zzdist
    fac[singular] = 1.0

    return fac.prod(axis=-1)

class BiasZ:
    '''  Define functions for b1(z) for different surveys  '''
    def __init__(self, cosmo):
        self.cosmo = cosmo

    def const(self,  b = 1.):
        return lambda z: b

    def linear(self, b=2., z0 = 1.):
        return lambda z: self.cosmo.Dgrow(z=z0)*b/self.cosmo.Dgrow(z=z)

    def lsst(self, b, z0 = 1):
        return lambda z: 1 + 0.84*z
    
class DnDz:
    '''  Define dn/dz functions for different surveys per steradian'''    

    def __init__(self):
        self.arcm = np.pi/180./60.
        self.ster = self.arcm ** -2.
    
    def lsst(self, z, ilim = 25.3, N0 = 46.):
        z0 = 0.0417 * ilim -0.744	# LSST science book (3.8)
        Ng = N0 * 10.0**(0.3*(ilim-25))	# LSST science book (3.7)
        zz = z/z0
        f = 0.5*zz**2*np.exp(-zz)/z0

        return Ng*f 

    def lsstfaint(self, z, ilim = 24.3, N0 = 46.):
        z0 = 0.0417 * ilim -0.744	# LSST science book (3.8)
        Ng = N0 * 10.0**(0.3*(ilim-25))	# LSST science book (3.7)
        zz = z/z0
        f = 0.5*zz**2*np.exp(-zz)/z0

        return Ng*f 
    
    def des(self, z, ilim= 23.50, N0 = 46.):
        z0 = 0.0417 * ilim -0.744	## Functional forms same as LSST
        Ng = N0 * 10.0**(0.3*(ilim-25))	
        zz = z/z0
        f = 0.5*zz**2*np.exp(-zz)/z0

        return Ng*f

    def desdeep(self, z, ilim= 24.50, N0 = 46.):
        z0 = 0.0417 * ilim -0.744	## Functional forms same as LSST
        Ng = N0 * 10.0**(0.3*(ilim-25))	
        zz = z/z0
        f = 0.5*zz**2*np.exp(-zz)/z0

        return Ng*f

    def subaru(self, z, ilim= 25.90, N0 = 46.):
        z0 = 0.0417 * ilim -0.744	## Functional forms same as LSST
        Ng = N0 * 10.0**(0.3*(ilim-25))	
        zz = z/z0
        f = 0.5*zz**2*np.exp(-zz)/z0

        return Ng*f 

    def func(self, survey = 'lsst'):
        if survey == 'lsst':
            return lambda z: self.lsst(z)
        if survey == 'lsstfaint':
            return lambda z: self.lsstfaint(z)
        elif survey == 'des':
            return lambda z: self.des(z)
        elif survey == 'desdeep':
            return lambda z: self.desdeep(z)
        elif survey == 'subaru':
            return lambda z: self.subaru(z)

def nbar(zmin, zmax, dndz=None, survey='lsst'):
    """   
    Returns nbar, in galaxies per arcmin^2, for a given survey or dndz function within zmin<z<zmax, using DnDz class
    """
    zz = np.linspace(zmin, zmax, 500)

    if dndz is None:
        dndz = DnDz()

        if survey == "lsst":
            ff = np.trapz(y = dndz.lsst(zz), x = zz)
       
        elif survey == "lsstfaint":
            ff = np.trapz(y = dndz.lsstfaint(zz), x = zz)
        
        elif survey == "des":
            ff = np.trapz(y = dndz.des(zz), x = zz)
        
        elif survey == "desdeep":
            ff = np.trapz(y = dndz.desdeep(zz), x = zz)
        
        elif survey == "subaru":
            ff = np.trapz(y = dndz.subaru(zz), x = zz)
        
        else:
            raise ValueError("%s is unidentified. \n Specify survey from lsst, des or subaru. Or give your own dndz" % survey)
            
    else:
        ff = np.trapz(y = dndz(zz), x = zz)

    return ff
    
def parse_arguments(Files = None, CosmoPar = None):
   ''' Parse the command line arguments ''' 
   
   if CosmoPar is None:
       raise ValueError("Provide cosmological params to parse_arguments.")

   parser = argparse.ArgumentParser()

   parser.add_argument('--beam',      default = 1,            type = float)
   parser.add_argument('--temp',      default = 1,            type = float)
   parser.add_argument('--zg',        default = 1,            type = float,  help ="Redshift of source")
   parser.add_argument('--dz',        default = 0.5,          type = float,  help ="Redshift bin around the source")
   parser.add_argument('--Nsource',   default = 50,           type = float,  help ="Number of sources per arcmin")
   parser.add_argument('--b1',        default = None,         type = float,  help ="Bias of the sources")
   parser.add_argument('--M',         default = CosmoPar.M,   type = float,  help ="Matter density")
   parser.add_argument('--L',         default = CosmoPar.L,   type = float,  help ="Dark energy")
   parser.add_argument('--H0',        default = CosmoPar.H0,  type = float,  help ="Hubble")
   parser.add_argument('--pkfile',    default = Files.pkfile,                help ="File with linear power spectrum")
   parser.add_argument('--outpath',   default = Files.outpath,               help ="Where to save output")
   parser.add_argument('--classfile', default = Files.classfile,             help ="cl output from class")
   parser.add_argument('--errorfile', default = None,                        help ="spectra file with errors and noises")
   parser.add_argument('--noisefile', default = None,                        help ="File with instrumental and Fisher Noise")
   parser.add_argument('--survey',    default = 'lsst',                      help ="Which survey to use dndz from")

   args        = parser.parse_args()
   args.cosmo  = cosmology.Cosmology(M = args.M, L = args.L, pfile = args.pkfile)

   if args.b1 is None:
       args.b1 = 2.*args.cosmo.Dgrow(z = 1.)/args.cosmo.Dgrow(z = args.zg)

   ## Define file formats here
   if args.noisefile is None:
       try:
           args.noisefile = Files.noisefile%(args.beam*10, args.temp*10)
       except:
           args.noisefile = Files.noisefile
   
   if args.errorfile is None:
       try:
           args.errorfile = Files.errorfile%(args.zg*100, args.b1*10, args.beam*10, args.temp*10)
       except:
           args.errorfile = Files.errorfile

   ## Print them for runtime check
   print("beam = ",  10.*args.beam)
   print("dT = ",    10.*args.temp)
   print("zg = ",        args.zg)
   print("b1 = ",        args.b1)
   print("Noisefile = ", args.noisefile)
   print("Errorfile = ", args.errorfile)

   return args

