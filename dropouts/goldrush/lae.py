import  os
import  matplotlib         as      mpl  
import  numpy              as      np
import  pandas             as      pd
import  matplotlib.pyplot  as      plt 

from    utils              import  pprint
from    cosmo              import  cosmo
from    params             import  get_params
from    scipy.interpolate  import  interp1d
from    utils              import  latexify
from    contamination      import  get_contamination
from    cosmo              import  cosmo


params = get_params()

def get_betas(dropband, depth, printit=False):
    '''    
    NOTE:  No zy for Ultra Deep. 
    '''
    import  os
    import  glob

    from    specs            import  samplestats
    from    selection_box    import  detection_bands


    root          = os.environ['LBGCMB']
    
    ##  Loads from (hsc) specs.py
    specs         = samplestats()

    ##  Catalogues labelled by bands used in selection criteria. 
    band2cat      = {'g': 'gri', 'r': 'riz', 'i': 'izy', 'z': 'zY'}

    if band2cat[dropband]  == "zy":     
        columns   =  pd.read_csv(root + "/dropouts/goldrush/cats/column_zy.cat", header=None, names=["colname"], delim_whitespace=True)
        appradii  =  '10'

    else:
        columns   =  pd.read_csv(root + "/dropouts/goldrush/cats/column.cat", header=None, names=["colname"], delim_whitespace=True)
        appradii  =  '15'

    
    ##  Define magnitude type to be retrieved. 
    mags          = [x + appradii     for x in ['gmag_aperture', 'rmag_aperture', 'imag_aperture', 'zmag_aperture', 'ymag_aperture']]
    ##  mags      = [x + 'cmodel_mag' for x in ['g', 'r', 'i', 'z', 'y']]
    
    ##  Columns to be retrieved from GoldRush catalogues. 
    names         = ['object_id', 'ra', 'dec'] + mags

    ##  Get all the fields for given depth.  
    files         = glob.glob(root + "/dropouts/goldrush/cats/%s/%s/*" % (depth, band2cat[dropband]))

    counts        =  []
    total_area    = 0.0

    ##  Loop over fields. 
    for file in files:                                                                                                              
        mid         = file.split(".")[-2]
        field       =  mid.split('_' + depth + '_')[-1]
        
        ##  Area [deg^2] of e.g. COSMOS 'D'
        area        = specs['areas'][depth][field]                                                                               
        total_area += area

        ##  Read and drop rows with undefined magnitudes. 
        input       = pd.read_csv(file, header=None, names=columns["colname"], delim_whitespace=True)
        input       = input.replace(99., np.nan)

        ##  Drop any with magnitude fault. 
        input       = input.dropna(axis=0, how='any', thresh=None, subset=[x + 'mag_aperture' + appradii for x in ['g', 'r', 'i', 'z', 'y']],\
                                   inplace=False)
                
        ##  Copy the necessary part of the dataframe. 
        data        = input[names].copy()

        pl.scatter(data['zmag_aperture' + appradii], data['imag_aperture' + appradii] - data['zmag_aperture' + appradii],\
                   c='k', marker='x', rasterized=True, s=7, alpha=0.4)

        ##  data['lsst-beta'] = 5.3 * (data['imag_aperture' + appradii] - data['zmag_aperture' + appradii]) - 2.04

    zs   = np.arange(20., 30., 0.1)
    
    beta = -2.0
    cs   = np.ones_like(zs) * (beta + 2.04) / 5.3

    pl.plot(zs, cs)

    pl.xlim(22.0, 26.0)
    pl.ylim(-1.0,  2.0)

    pl.xlabel(r'$z$')
    pl.ylabel(r'$(i-z)$')

    pl.show()
 
if __name__ == "__main__":
    import  collections
    import  pylab          as  pl
    import  astropy.units  as  u


    print('\n\nWelcome to the HSC LAE calculator.')
    
    for zee in [3.0, 4.0, 5.0]:
        DL               = cosmo.luminosity_distance(zee).value  ## [Mpc]  
        DL              *= 1.e5

        ms               = np.arange(21.5, 30.5, 0.5)
        Ms               = ms + 2.5 * np.log10(1. + zee) - 5.0 * np.log10(DL)

        xlya             = np.zeros_like(Ms)
  
        xlya[Ms < -20.3] = 1.09 + 0.047 * Ms[Ms < -20.3]
        xlya[Ms > -20.3] = 5.46 +  0.26 * Ms[Ms > -20.3]

        ##  pl.plot(Ms, ms)
        pl.plot(ms, xlya, label=r'$z=%.1lf$' % zee)

    pl.xlabel(r'$m_{5 \sigma}$')
    pl.ylabel(r'Dropout fraction with EW > 50 angstroms')

    pl.xlim(22.5, 27.)
    pl.ylim(0.0, 1.0)
    pl.legend()
    pl.show()
    
    ##  get_betas('g', 'D', printit=True)
    
    print("\n\nDone.\n\n")
