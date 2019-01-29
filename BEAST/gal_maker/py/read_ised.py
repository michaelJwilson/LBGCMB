import sys 
import math
import array, os, re
import numpy              as      np
import astropy.constants  as      const

from    astropy           import  units  as u


# some useful constants
c         = const.c.value       # speed of light (m/sec)
m_per_au  = 1.49598e11          # meters per astronomical unit
au_per_pc = 3600 * 180 / np.pi  # AUs per parsec

def to_hertz(to_convert, units='a'):
    '''
    res = ezgal.utils.to_hertz( to_convert, units='Angstroms' )
    Converts the given wavelength (in the given units) to hertz.
    :param to_convert: The wavelength to convert
    :param units: The units the wavelength is in
    :type to_convert: int, float
    :type units: string
    :returns: The converted frequency
    :rtype: float
    :Example:
    >>> import ezgal
    >>> ezgal.utils.to_hertz( 1000, units='a' )
    2997924580000000.0
    .. seealso::
    see :func:`ezgal.utils.to_meters` for list of available units
    Also see :func:`ezgal.utils.to_lambda`
    '''

    return (c / to_meters(1.0, units=units)) / np.asarray(to_convert)

def to_meters(to_convert, units='a'):
    '''
    res = ezgal.utils.to_meters( to_convert, units='a' )
    Converts a length from the given units to meters
    :param to_convert: The length to convert
    :param units: The units to convert the length to
    :type to_convert: int, float
    :type units: string
    :returns: The converted length
    :rtype: int, float
    :Example:
    >>> import ezgal
    >>> ezgal.utils.to_meters( 1e10, units='a' )
    1.0

    **units**
    Available units are (case insensitive):
    ================= ====================
           Name             Units
    ================= ====================
    a,angstroms       Angstroms
    nm,nanometers     Nanometers
    um,microns        Microns
    mm,millimeters    Millimeters
    cm,centimeters    Centimeters
    m,meters          Meters
    km,kilometers     Kilometers
    au                Astronomical Units
    pc,parsecs        Parsecs
    kpc, kiloparsecs  Kiloparsecs
    mpc, megaparsecs  Megaparsecs
    ================= ====================
    .. seealso:: :func:`ezgal.utils.convert_length`
    '''

    units = units.lower()
    to_convert = np.asarray(to_convert)

    if units == 'angstroms' or units == 'a': return to_convert / 1e10
    if units == 'nanometers' or units == 'nm': return to_convert / 1e9
    if units == 'microns' or units == 'um': return to_convert / 1e6
    if units == 'millimeters' or units == 'mm': return to_convert / 1e3
    if units == 'centimeters' or units == 'cm': return to_convert / 1e2
    if units == 'meters' or units == 'm': return to_convert
    if units == 'kilometers' or units == 'km': return to_convert * 1000.0
    if units == 'au': return to_convert * m_per_au
    if units == 'parsecs' or units == 'pc':
        return to_convert * m_per_au * au_per_pc
    if units == 'kilparsecs' or units == 'kpc':
        return to_convert * m_per_au * au_per_pc * 1000.0
    if units == 'megaparsecs' or units == 'mpc':
        return to_convert * m_per_au * au_per_pc * 1e6

    raise NameError('Units of %s are not recognized!' % units)

def convert_length(to_convert, incoming='m', outgoing='a'):
    '''
    res = ezgal.utils.convert_length( to_convert, incoming='m', outgoing='a' )
    converts a length from the incoming units to the outgoing units.
    :param to_convert: The length to convert
    :param incoming: The units to convert the length from
    :param outgoing: The units to convert the length to
    :type to_convert: int, float
    :type incoming: string
    :type outgoing: string
    :returns: The converted length
    :rtype: int, float
    :Example:
    >>> import ezgal
    >>> ezgal.utils.convert_length( 1, incoming='pc', outgoing='au' )
    206264.80624709636
    .. seealso:: see :func:`ezgal.utils.to_meters` for available units.
    '''

    return to_meters(to_convert, units=incoming) / to_meters(1.0,
                                                             units=outgoing)

def _read_binary(fhandle, type='i', number=1, swap=False):
    '''
    res = ezgal.utils._read_binary( fhandle, type='i', number=1, swap=False )
    reads 'number' binary characters of type 'type' from file handle 'fhandle'
    returns the value (for one character read) or a numpy array.
    set swap=True to byte swap the array after reading
    '''

    
    if (sys.version_info >= (3, 0)) & (type == 'c'):
       ##  unsigned char in python 2.  
       ##  https://docs.python.org/2/library/array.html
       ##  https://docs.python.org/3/library/array.html
       ##  type  =  'B'  ##  unsigned char in python 3.
       ##  type  =  'b'  ##  signed char in python 3.
    
       import warnings 

       type = 'B'

       warnings.warn('Reassigning unsigned char type (c to B) as per python 3.')


    arr = array.array(type)
    arr.fromfile(fhandle, number)

    if swap: 
        arr.byteswap()

    ## if sys.version_info >= (3, 0):
        ##  https://github.com/desihub/redrock/blob/master/py/redrock/utils.py                                                                                                                                                                                                     
        ##  if not array.dtype.isnative:                                                                                                                                                                                                                                             
        ##    return data.byteswap().newbyteorder()   

    if len(arr) == 1:
        return  arr[0]

    else:
        return np.asarray(arr)

def read_ised(file):
    '''
    Sourced from:  https://github.com/cmancone/easyGalaxy/blob/master/ezgal/utils.py
 
    ( seds, ages, vs ) = ezgal.utils.read_ised( file )
    Read a bruzual and charlot binary ised file.
    :param file: The name of the ised file
    :type file: string
    :returns: A tuple containing model data
    :rtype: tuple
    .. note::
    All returned variables are numpy arrays.  ages and vs are one dimensional arrays, and seds has a shape of (vs.size,ages.size)
    **units**
    Returns units of:
    =============== ===============
    Return Variable   Units
    =============== ===============
    seds            Ergs/s/cm**2/Hz
    ages            Years
    vs              Hz
    =============== ===============
    '''

    if not (os.path.isfile(file)):
        raise ValueError('The specified model file was not found!')

    print('Reading .ised:  %s' % str(file))

    # open the ised file
    fh    =  open(file, 'rb')

    # start reading
    junk  = _read_binary(fh)
    nages = _read_binary(fh)

    print(junk)

    # first consistency check
    if nages < 1 or nages > 2000:
        raise ValueError(
            'Problem reading ised file - unexpected data found for the number of ages!')

    # read ages
    ages = np.asarray(_read_binary(fh, type='f', number=nages))

    # read in a bunch of stuff that I'm not interested in but which I read like this to make sure I get to the right spot in the file
    junk = _read_binary(fh, number=2)
    iseg = _read_binary(fh, number=1)

    if iseg > 0: 
        junk = _read_binary(fh, type='f', number=6 * iseg)

    junk = _read_binary(fh, type='f', number=3)
    junk = _read_binary(fh)
    junk = _read_binary(fh, type='f')
    junk = _read_binary(fh, type='c', number=80)
    junk = _read_binary(fh, type='f', number=4)
    junk = _read_binary(fh, type='c', number=160)
    junk = _read_binary(fh)
    junk = _read_binary(fh, number=3)

    # read in the wavelength data
    nvs = _read_binary(fh)

    # consistency check
    if nvs < 10 or nvs > 12000:
        raise ValueError('Problem reading ised file - unexpected data found for the number of wavelengths!')

    # read wavelengths and convert to frequency (comes in as Angstroms)
    # also reverse the array so it will be sorted after converting to frequency
    ls   = _read_binary(fh, type='f', number=nvs)[::-1]

    # create an array for storing SED info
    seds = np.zeros((nvs, nages))

    # now loop through and read in all the ages
    for i in range(nages):
        junk = _read_binary(fh, number=2)
        nv   = _read_binary(fh)

        if nv != nvs:
            raise ValueError(
                'Problem reading ised file - unexpected data found while reading seds!')

        seds[:, i] = _read_binary(fh, type='f', number=nvs)[::-1]

        nx   = _read_binary(fh)
        junk = _read_binary(fh, type='f', number=nx)

    # now convert the seds from Lo/A to ergs/s/Hz
    seds   *= 3.826e33
    seds   *= ls.reshape((nvs, 1))**2.0
    seds   /= convert_length(c, outgoing='a') 

    # convert from ergs/s/Hz to ergs/s/Hz/cm^2.0 @ 10pc
    # seds /= 4.0 * np.pi * convert_length(10, incoming='pc', outgoing='cm')**2.0

    ##  Convert ages to Gyr                                                                                                     
    ages   /= 1.e9

    vs      = to_hertz(ls)

    fh.close()

    # sort in frequency space
    sinds = vs.argsort()

    return  (seds[sinds,:], ages, vs[sinds], ls[sinds])


if __name__ == '__main__':
  import os 
  import glob
  import pylab as pl


  root   =  os.environ['BEAST']
  
  ##  fname  = '/global/homes/m/mjwilson/desi/BEAST/gal_maker/dat/BC03/bc2003_hr_m62_salp_ssp.ised'
  ##  fname  = '/global/homes/m/mjwilson/desi/BEAST/gal_maker/dat/BC03/bc2003_hr_m22_chab_ssp.ised'
  ##  fname  = '/global/homes/m/mjwilson/desi/BEAST/gal_maker/dat/BOWLER/ised/bc2003_lr_m42_chab_tau_0_05.ised'

  files      = glob.glob(root + '/gal_maker/dat/BOWLER/ised/*')
    
  for fname in files:
    title = fname.split('/')[-1]

    ##  ergs/s/Hz, Giga Years, Hertz, Angstroms. 
    seds, ages, vs, ls = read_ised(fname)

    for ii, sed in enumerate(seds.T):
     if (ages[ii] > 0.01) & (ii % 25 == 0):    
        ##  vs * Fv = ls * Fl.  Fig. 13.5 of Cosmological Physics.
        pl.loglog(ls, vs * sed, label=ages[ii])

    pl.xlim(3.e2,  3.e4)
    pl.ylim(1.e28, 1.e38)
    pl.title(title)

    pl.legend()
    pl.show()
    

