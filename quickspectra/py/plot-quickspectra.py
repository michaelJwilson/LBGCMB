from    __future__            import division

import  astropy.io.fits       as pyfits
import  pylab                 as pl
import  numpy                 as np
import  matplotlib.pyplot     as plt
import  os 

from    matplotlib.offsetbox  import TextArea, AnnotationBbox


def get_input(ax, redshift = 5.5):
    fname = os.environ['SPECTRA'] + '/BC03/' + 'spec-BC03-z%.1lf.dat' % redshift

    input = np.loadtxt(fname)
    
    ax.plot(input[:,0], np.tanh(input[:,1]), 'c', label='Input flux.')

    return

def get_quickspectra(redshift = 5.5, exposure = 9000):
    ##  spec-BC03-z4.5_exp900.fits
    fname                = './dat/BC03/spec-BC03-z%.1lf_exp%.0lf.fits' % (redshift, exposure)
    exposure            /= 60.*60.  ## Hours
    
    title                = fname.split('/')[-1]

    h                    = pyfits.open(fname)

    h.info()
    
    '''
    Filename: spectra.fits
    No.    Name         Type      Cards   Dimensions   Format
    0   PRIMARY     PrimaryHDU      36   ()      
    1   FIBERMAP    BinTableHDU     69   10R x 28C   [10A, 20A, 8A, K, K, K, K, 5E, 50A, J, J, J, J, J, J, E, D, D, D, D, D, D, D, D, E, E, J, J]   
    2   B_WAVELENGTH  ImageHDU        10   (2380,)   float32   
    3   B_FLUX      ImageHDU        11   (2380, 10)   float32   
    4   B_IVAR      ImageHDU        10   (2380, 10)   float32   
    5   B_MASK      CompImageHDU     12   (2380, 10)   int32   
    6   B_RESOLUTION  ImageHDU        11   (2380, 9, 10)   float32   
    7   R_WAVELENGTH  ImageHDU        10   (2116,)   float32   
    8   R_FLUX      ImageHDU        11   (2116, 10)   float32   
    9   R_IVAR      ImageHDU        10   (2116, 10)   float32   
    10  R_MASK      CompImageHDU     12   (2116, 10)   int32   
    11  R_RESOLUTION  ImageHDU        11   (2116, 9, 10)   float32   
    12  Z_WAVELENGTH  ImageHDU        10   (2399,)   float32   
    13  Z_FLUX      ImageHDU        11   (2399, 10)   float32   
    14  Z_IVAR      ImageHDU        10   (2399, 10)   float32   
    15  Z_MASK      CompImageHDU     12   (2399, 10)   int32   
    16  Z_RESOLUTION  ImageHDU        11   (2399, 11, 10)   float32  
    '''

    fig, ax = plt.subplots()
    
    pl.plot(h["B_WAVELENGTH"].data, np.tanh(h["B_FLUX"].data[0]), 'b')
    pl.plot(h["R_WAVELENGTH"].data, np.tanh(h["R_FLUX"].data[0]), 'r')
    pl.plot(h["Z_WAVELENGTH"].data, np.tanh(h["Z_FLUX"].data[0]), 'g')

    get_input(ax)

    try:
      get_redshift(ax, exposure)

      pl.title(r"DESI %s Spec. @ $z = %.2lf$" % (type.upper(), np.float(z[1:])))

    except:
        pl.title(r"DESI %s Spectra" % title)

    pl.xlim(3.e3, 1.e4)
    ## pl.xlim(h["B_WAVELENGTH"].data.min(), h["Z_WAVELENGTH"].data.max())

    pl.xlabel(r"$\lambda$ [$\AA$]")
    pl.ylabel(r"$\tanh(F_\lambda$ / $10^{-17}/$ ergs/$s$/cm$^2$/$\AA$)")

    pl.savefig('plots/quickspectra.pdf')

def get_redshift(ax, exposure):
    h = pyfits.open('../redrock/' + os.environ['FILE'] + '-zbest.fits')

    h.info()
    
    '''
    Filename: desispec-lrg-z0.7-zmag20.00_zbest.fits
    No.    Name         Type      Cards   Dimensions   Format
    0  PRIMARY     PrimaryHDU       4   ()      
    1  ZBEST       BinTableHDU     29   1R x 10C   [D, 10D, D, D, K, 6A, 1A, K, D, 8A]
    '''
    '''
    XTENSION= 'BINTABLE'           / binary table extension                         
    BITPIX  =                    8 / array data type                                
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                  143 / length of dimension 1                          
    NAXIS2  =                    1 / length of dimension 2                          
    PCOUNT  =                    0 / number of group parameters                     
    GCOUNT  =                    1 / number of groups                               
    TFIELDS =                   10 / number of table fields                         
    TTYPE1  = 'CHI2    '                                                            
    TFORM1  = 'D       '                                                            
    TTYPE2  = 'COEFF   '                                                            
    TFORM2  = '10D     '                                                            
    TTYPE3  = 'Z       '                                                            
    TFORM3  = 'D       '                                                            
    TTYPE4  = 'ZERR    '                                                            
    TFORM4  = 'D       '                                                            
    TTYPE5  = 'ZWARN   '                                                            
    TFORM5  = 'K       '                                                            
    TTYPE6  = 'SPECTYPE'                                                            
    TFORM6  = '6A      '                                                            
    TTYPE7  = 'SUBTYPE '                                                            
    TFORM7  = '1A      '                                                            
    TTYPE8  = 'TARGETID'                                                            
    TFORM8  = 'K       '                                                            
    TTYPE9  = 'DELTACHI2'                                                           
    TFORM9  = 'D       '                                                            
    TTYPE10 = 'BRICKNAME'                                                           
    TFORM10 = '8A      '                                                            
    EXTNAME = 'ZBEST   '  
    '''
    
    print
    print "Best fit    z:  %.3lf" % h[1].data["Z"][0] 
    print "Error in    z:  %.5le" % h[1].data["ZERR"][0]
    print
    print "Best fit  X^2:  %.3lf" % h[1].data["CHI2"][0]
    print "Best fit dX^2:  %.3lf" % h[1].data["DELTACHI2"][0]

    offsetbox = TextArea("Best fit z:  %.3lf +- %.1le\nExposure time: %.1lf hrs" % (h[1].data["Z"][0],\
                          h[1].data["ZERR"][0], exposure), minimumdescent=False)

    xy        = (5500., 1.5)
    xy        = (0.70, 0.82)
    
    ab        = AnnotationBbox(offsetbox, xy,
                               xybox=(0, 0),
                               xycoords='figure fraction',
                               boxcoords="offset points")
    
    ax.add_artist(ab)


if __name__ == "__main__":
    get_quickspectra()
