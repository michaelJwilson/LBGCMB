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
from    effective_depth    import  effective_depth
from    galaxy_frac        import  galaxy_frac
from    contamination      import  get_contamination


params = get_params()

def get_nbar_nocontam(band, depth='W', printit=False):
    '''                                                                                                                                                      
    Get the number counts per magnitude bin, and the contamination corrected counts.                                                                        
    '''

    from  goldrush.specs  import  samplestats


    magbins, pnbar, counts    = get_nbarbymag(band, depth, printit=False)

    stats                     = samplestats(printit=printit)
    zee                       = np.ceil(stats[band]['z'])     ## Round to the nearest integer.                                                               

    stats[band]['nbar_noint'] = sum(counts * (1. - get_contamination(magbins[:-1], zee=zee, depth=depth))) / stats['Total area'][depth]

    if printit:
      print  sum(stats[band]['counts'].values()), counts.sum(), sum(counts * (1. - get_contamination(magbins[:-1], zee=zee, depth=depth)))

    return  stats

def get_nbarbymag(dropband, depth, printit=False):
    '''    
    For selection bands, e.g. gri, and depth, [UD, D, W], 
    calculate the total number of objects with a magnitude less 
    than a given limit and the area of that field type. 
    
    From this, calculate the angular number density.  NOTE:  No zy for Ultra Deep. 
    '''
    import  os
    import  glob

    from    des              import  des_depths
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
    mags          = [x + appradii for x in ['gmag_aperture', 'rmag_aperture', 'imag_aperture', 'zmag_aperture', 'ymag_aperture']]
    
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
        input       = input.dropna(axis=0, how='any', thresh=None, subset=[detection_bands[dropband] + 'mag_aperture' + appradii], inplace=False)
                
        ##  Copy the necessary part of the dataframe. 
        data        = input[names].copy()

        if printit:
            print("\n\nHSC %s field with area %.6lf deg2 (from file %s). \n" % (field, area, file))                                                         

            pprint(data)                                                                                                                                 

            print("\n\nMagnitude limits of identified objects:  %.6lf < mag < %.6lf" % \
                 (data[detection_bands[dropband] + 'mag_aperture' + appradii].min(), data[detection_bands[dropband] + 'mag_aperture' + appradii].max()))     

        if input.empty == False:
           ##  Bands in mag. for the detection band, to calculate Ng( < mag). 
           magbins      = np.arange(18., 26.5, 0.01)                 

           ##  Catch all bin for objects brighter than 22nd mag and fainter than 26.5.
           magbins      = np.concatenate([np.array([0.]), magbins, np.array([30.])])
           
           ##  Defines the bin edges, including the rightmost edge, allowing for non-uniform bin widths.
           counts.append(np.histogram(data[detection_bands[dropband] + 'mag_aperture' + appradii], bins=magbins)[0])

    counts = np.array(counts)

    ##  Sum the binned histogram counts over all fields. 
    counts = counts.sum(axis=0)

    if not (depth is 'UD' and band2cat[dropband] is 'zy'):
      cumulative  = np.cumsum(counts)        ##  Cumulative counts of all bins up to mag. limit;    
      pnbar       = cumulative / total_area  ##  Projected nbar.  

      ##  Check.                                                                                                                                                                                                                  
      print("\n\nTotal  counts @ depth %s (%.2lf):  %.3lf \t %.3lf"  % (depth, effective_depth(dropband=dropband, depth=depth), cumulative[-1],\
                                                                        specs[dropband]['total-%scounts' % {'W': '', 'D': 'deep', 'UD': 'udeep'}[depth]]))
      print("Total     area:  %.3lf sq. deg. (%.3lf)"                % (total_area,     specs['Total area'][depth]))
      print("Effective nbar:  %.3lf per sq. deg."                    % (cumulative[-1] / total_area))

      return  magbins, pnbar, counts

    else:
      raise  ValueError("\n\nThe combination of UD and zy band selection is not available.\n\n")
        
def plot_ilims(results, plot_des = False, plot_hsc = True):
    import  os
    import  pylab              as      pl
    import  matplotlib.pyplot  as      plt
    
    from    des                import  des_depths
    from    matplotlib.pyplot  import  text
    from    specs              import  samplestats
    from    selection_box      import  detection_bands


    latexify(fig_width=None, fig_height=None, columns=2, equal=False)
    
    colors = ['b', 'r', 'indigo']

    ##  Load Goldrush basic stats. 
    stats  =  samplestats(printit=False)
    
    ##  Create figure.                                                                                                                                       
    fig, axarray    = plt.subplots(1, 1, sharey=False)
    fig.set_size_inches(6.5, 3.5)

    ##  Catch a one plot call. 
    if not isinstance(axarray, np.ndarray):
      axarray = np.array([axarray])

    for k, dropband in enumerate(results.keys()):        
        for ldepth, label in zip(['D'], ['GOLDRUSH Deep']):
            magbins, pnbar, counts = results[dropband][ldepth]

            ##  Plot ang_nbar against limiting mag. (rightmost edge).                                                                                                                                         
            axarray[k].semilogy(magbins[1:], pnbar, '-',  label=r'$g-$' + label, markersize=3, lw=1.)

        for contamination, color in zip(['D', 'W'], ['dodgerblue', 'indigo']):
          rate = get_contamination(magbins[1:], round(stats[dropband]['z']), contamination)

          axarray[k].semilogy(magbins[1:], pnbar * (1. - rate), '--', label='Less %s' % contamination + r' ($\simeq$' + '%.2lf) interlopers' % effective_depth(dropband, contamination),\
                                                                      markersize=3, lw=1., c=color, dashes=[3, 1], alpha=0.4)

        if plot_des:
          ##  DES depths/yr estimate                                                                                                                    
          depths = des_depths()

          ##  Plot the magnitude limit of DES SV. 
          axarray[k].axvline(depths['SV'][detection_bands[dropband]], c='k', linestyle='-', label='DES SV', lw = 0.5)
        
          ##  and the magnitude limits per year (estimated for year greater than one).
          for year in np.arange(1, 6, 1):        
            axarray[k].axvline(depths['Y' + str(year)][detection_bands[dropband]], c='k', linestyle='-', label='', lw = 0.5)

        root = os.environ['LBGCMB']
        data = np.loadtxt(root + "/dropouts/nz/schechter/dat/schechter_estimate_%s_dropouts.txt" % dropband)

        ##  0.59 completeness frac.
        axarray[k].semilogy(data[:,0], data[:,1], 'k', label='Best-fit UV Schechter fn.', alpha=0.5)

        ##  0.59 completeness frac.  Galaxies only. 
        axarray[k].semilogy(data[:,0], galaxy_frac(data[:,0]) * data[:,1], 'k--', label='Best-fit galaxy UV Schechter fn.', dashes=[3,1])
        
        ymax = pnbar[magbins < 24.5].max()
        
        axarray[k].set_xlim(23.0,       26.5)  ## [magbins[1], magbins[-2]] 
        axarray[k].set_ylim([10., 8.e3])

        axarray[k].set_xlabel(r'$%s_{\rm{AB}}$'     %  detection_bands[dropband])
        axarray[k].set_ylabel(r'$N(<%s)$ / deg$^2$' %  detection_bands[dropband])

        ## axarray[k].get_yaxis().get_major_formatter().set_scientific(False)
        ## axarray[k].ticklabel_format(style='sci', scilimits=(-1, 2))

    pl.legend(loc=4)
    pl.savefig('plots/hsc_icut.pdf', bbox_inches='tight')


if __name__ == "__main__":
    import  collections
    import  pylab as pl

    print('\n\nWelcome to the HSC dropout ilim calculator.')

    
    magbins, pnbar, counts = get_nbarbymag('r', 'D', printit=False)

    pl.semilogy(magbins[:-1], pnbar)
    pl.xlim(20., 27.0)
    pl.show()
    
    stats  = get_nbar_nocontam('r', depth='W', printit=False)                                                                                           
    pprint(stats) 

    results = collections.OrderedDict()

    for k, dropband in enumerate(['g']):
      ##  'UD':  get_nbarbymag(dropband, 'UD', printit=False)
      results[dropband] = {'W': get_nbarbymag(dropband, 'W', printit=False), 'D':  get_nbarbymag(dropband, 'D', printit=False)}

    print('')

    plot_ilims(results)
    
    print("\n\nDone.\n\n")
