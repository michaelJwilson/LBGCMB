import  numpy              as np
import  pylab              as pl
import  pandas             as pd
import  matplotlib.cm      as cm
import  matplotlib.pyplot  as plt
import  scipy.interpolate  as interpolate
import  numpy.ma           as ma

'''
-- GOLDRUSH selection criteria --
Detection band(s):  We select sources with signal-to-noise ratio (S/N) > 5 within 1.5'' diameter apertures:                                                
g-dropouts:  i                                                                                                                                    
r-dropouts:  z                                                                                                                                   
i-dropouts:  z (and 4.0 sigma in y)                                                                                                                 
z-dropouts:  y                                                                                                                                    

To remove low-z source contamination, we also require that sources be undetected (< 2 sigma) within 1.5'' diameter apertures:                            
r-dropouts: in g                                                                                                                           
i-dropouts: in g, r                                                                                                                              
z-dropouts: in g, r, i                                                                                                                          

Returns true if galaxy meets the 5 sig limit; 
'''
'''
## Assumes best 5 sigma depth (UD and best of all fields) for each band:                                      
g5s = 27.15
r5s = 26.84
i5s = 26.53
z5s = 26.10
y5s = 25.28
'''
##  Shallower W survey:
g5s = 26.40
r5s = 25.90
i5s = 25.70
z5s = 25.00
y5s = 24.25

## Band for (5 sig.) detection; set to HSC Goldrush selection.
detection_bands  = {'g':  'i', 'r':  'z', 'i':  'z', 'z':  'y'}

def Subaru_significancecut(mags, dropband='g', sigcut = True, check_dropbands = True, masked_array = False):
    ''' 
    Returns (zs, mass, Luv, u, g, r, i, z, y) as an array from the input 
    spliced by selection criteria of dropband or masked array.
    '''
    zs   = mags['zee']
    mass = mags['mass']
    Luv  = mags['Luv']

    u, g, r, i, z, y = mags['u'], mags['g'], mags['r'], mags['i'], mags['z'], mags['y']

    ## Given dropband, assign the selection criteria to crit. 
    if dropband == 'g':
        ## 5 sigma detection in i for g-drops. 
        crit  = i < i5s        
                                                                                                                         
    elif dropband == 'r':
        ## 5 sigma detection in z for r-drops.
        crit  = z < z5s                                                                                                                 
        
        if check_dropbands == True:
            ## Remove low-z source contamination by ensuring sources undetected (< 2 sigma) in g.
            crit &= g > (g5s - 2.5 * np.log10(2./5.)) 

    elif dropband == 'i':
        ## 5 sigma detection in z for i-drops.
        crit  = z <  z5s

        ## And 4.0 sigma detection in y. 
        crit &= y < (y5s - 2.5 * np.log10(4./5.))

        if check_dropbands == True:
            ## Remove low-z source contamination by ensuring sources undetected (< 2 sigma) in g+r
            crit &= (g > (g5s - 2.5 * np.log10(2./5.))) & (r > (r5s - 2.5 * np.log10(2./5.)))                                             

    elif dropband == 'z':
        ## 5 sigma detection in y for z-drops.
        crit = y < y5s

        if check_dropbands == True:
            ## Remove low-z source contamination by ensuring sources undetected (< 2 sigma) in g+r+i
            crit &= (g > (g5s - 2.5 * np.log10(2./5.))) & (r > (r5s - 2.5 * np.log10(2./5.))) & (i > (i5s - 2.5 * np.log10(2./5.)))

    else:
        raise ValueError("Available HSC dropbands are g, r, i, z.")

    if masked_array:
        if sigcut:
            ## This doesn't look right. 
            return crit
        
        else:
            return ma.make_mask(np.ones_like(zs))
            
    else:
        if sigcut:
            return (zs[crit], mass[crit], Luv[crit], u[crit], g[crit], r[crit], i[crit], z[crit], y[crit])

        else:
            return (zs, mass, Luv, u, g, r, i, z, y)

def Subaru_colourcut(mags, band='g', nocolourcut = False, mask = False):
    ''' 
    Returns zs and Luvs of input objects whose magnitudes satisfy GOLDRUSH colour cuts 
    '''                                                    

    zs   = mags['zee']
    mass = mags['mass']
    Luv  = mags['Luv']

    u    = mags['u']
    g    = mags['g']
    r    = mags['r']
    i    = mags['i']
    z    = mags['z']
    y    = mags['y']

    if band == 'g':
        crit  = ((g - r) > 1.0) & ((r - i) < 1.0) & ((g - r) > 1.5*(r - i) + 0.8) 
        
    elif band == 'r':
        crit  = ((r - i) > 1.2) & ((i - z) < 0.7) & ((r - i) > 1.5*(i - z) + 1.0)

    elif band == 'i':
        crit  = ((i - z) > 1.5) & ((z - y) < 0.5) & ((i - z) > 2.0*(z - y) + 1.1)

    elif band == 'z':
        crit  = (z - y > 1.6)

    else:
        raise ValueError("Available bands are g, r, i, z.")

    if mask == True:
        return crit

    else:
        if nocolourcut == True:
            return (zs, mass, Luv, u, g, r, i, z, y)
 
        else:
            return (zs[crit], mass[crit], Luv[crit], u[crit], g[crit], r[crit], i[crit], z[crit], y[crit])

def get_SubaruSelectionfn(printit = True, plotit = True):
    '''
    Based on Bruzual & Charlot SEDS (Le Phare output of .phys and .output) and Subaru colour cut criteria, return a (selection) function specified 
    by redshift and UV absolute magnitude arguments. 
    '''

    EBVs = get_pEBV()
    Hs   = []

    results = {'g': [], 'r': [], 'i': [], 'z': []}

    for EBV in np.array([-0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45]):
       mags          = pd.read_csv(os.environ['CMBLBG'] + 'dropouts/test_ebv_%.3lf.dat' % EBV,\
                            names=["zee", "age", "mass", "Luv", "u", "g", "r", "i", "z", "y"], delim_whitespace=True)  ## stellar mass. 
    
       ## Generate result object.
       result        = np.zeros_like(mags['z'])

       ## Allow for arbitrary flux scaling.                                                                                                                
       flux_scaling  = 1.0 
    
       mags['mass'] *= flux_scaling
       mags['Luv']  *= flux_scaling

       mags['u']    -= 2.5 * np.log10(flux_scaling)
       mags['g']    -= 2.5 * np.log10(flux_scaling)
       mags['r']    -= 2.5 * np.log10(flux_scaling)
       mags['i']    -= 2.5 * np.log10(flux_scaling)
       mags['z']    -= 2.5 * np.log10(flux_scaling)
       mags['y']    -= 2.5 * np.log10(flux_scaling)
       
       mags['g_sigcut']   = Subaru_significancecut(mags, band='g',      sigcut =  True, check_dropbands = True, mask = True)
       mags['g_colorcut'] =       Subaru_colourcut(mags, band='g', nocolourcut = False, mask = True)
    
       mags['r_sigcut']   = Subaru_significancecut(mags, band='r',      sigcut =  True, check_dropbands = True, mask = True)
       mags['r_colorcut'] =       Subaru_colourcut(mags, band='r', nocolourcut = False, mask = True)
       
       mags['i_sigcut']   = Subaru_significancecut(mags, band='i',      sigcut =  True, check_dropbands = True, mask = True)
       mags['i_colorcut'] =       Subaru_colourcut(mags, band='i', nocolourcut = False, mask = True)

       mags['z_sigcut']   = Subaru_significancecut(mags, band='z',      sigcut =  True, check_dropbands = True, mask = True)
       mags['z_colorcut'] =       Subaru_colourcut(mags, band='z', nocolourcut = False, mask = True)

       if printit == True:
           print "\n\nObserved magnitudes layout (flux scaling of %.3lf)" % flux_scaling
           print mags

       flux_scalings     = np.logspace(0.0, 6.0, 1000)
       '''
       ## Bin selection in UV luminosity and redshift. 
       H, xedges, yedges = np.histogram2d(np.array([mags['zee'].values for x in flux_scalings]).flatten(),\
                                          np.array([np.log10(x * mags['mass'].values) for x in flux_scalings]).flatten(),\
                                          bins=[30, 20])
       '''
       H, xedges, yedges = np.histogram2d(np.array([mags['zee'].values for x in flux_scalings]).flatten(),\
                                          np.array([np.log10(x * mags['Luv'].values) for x in flux_scalings]).flatten(),\
                                          bins=[30, 20])
       
       Hs.append(H)

       if printit == True:
           print "\n\n%.3lf < z < %.2lf: "           % (xedges.min(), xedges.max())
           ## print "%.3lf < log10(Mass) < %.3lf: "  % (yedges.min(), yedges.max()), "\n\n"
           print "%.3lf < log10(Luv) < %.3lf: "      % (yedges.min(), yedges.max()), "\n\n"
    
       for i, x in enumerate(results.keys()):
         HH = np.zeros_like(Hs[0])  

         for scaling in flux_scalings:
            zs                            = mags['zee'].values

            mass                          = np.log10(mags['mass'].values * scaling)
            Luv                           = np.log10(mags['Luv'].values  * scaling)

            u                             = mags['u'].values - 2.5 * np.log10(scaling)
            g                             = mags['g'].values - 2.5 * np.log10(scaling)
            r                             = mags['r'].values - 2.5 * np.log10(scaling)
            i                             = mags['i'].values - 2.5 * np.log10(scaling)
            z                             = mags['z'].values - 2.5 * np.log10(scaling)
            y                             = mags['y'].values - 2.5 * np.log10(scaling)

            (zs, mass, Luv, u, g, r, i, z, y) = Subaru_significancecut({'zee': zs, 'mass': mass, 'Luv': Luv, 'u': u, 'g': g, 'r': r, 'i': i, 'z': z, 'y':y},\
                                                                         band=x, sigcut = True, check_dropbands = True, mask = False)
            
            ## Apply subaru colour cuts to catalogue, given their apparent magnitudes. 
            (zs, mass, Luv, u, g, r, i, z, y) = Subaru_colourcut({'zee': zs, 'mass': mass, 'Luv': Luv, 'u': u, 'g': g, 'r': r, 'i': i, 'z': z, 'y':y},\
                                                                   band=x, nocolourcut = False, mask = False)

            ## Plot selection as 2D histogram in z and L_uv.     
            ## TH, xxedges, yyedges        = np.histogram2d(zs, mass, bins=[xedges, yedges])
            TH, xxedges, yyedges           = np.histogram2d(zs, Luv, bins=[xedges, yedges])

            HH += TH
            
         results[x].append(HH)    

    H  = np.zeros_like(Hs[0])

    for x in Hs:
        H += EBVs['p(z=3)'] * x

    HH = np.zeros_like(Hs[0]) 

    for x in results.keys():
        x

    ## Selection: (Models retained / Models input).                                                                                                         
    HH /= H.astype(float)
        
    if plotit == True:
        pl.clf()

        plt.imshow(HH.T, origin='lower', extent=(xedges.min(), xedges.max(), yedges.min(), yedges.max()), interpolation='bicubic', aspect='auto')
            
        plt.xlabel(r'$z$')

        ## plt.ylabel(r'log$_{10}$(Stellar Mass)')
        plt.ylabel(r'log$_{10}$(L_uv)')

        ## ax = plt.gca()
        ## ax.set_aspect(0.1)
                
        plt.colorbar()

        ## pl.savefig('selection_%s.pdf' % x)
        pl.savefig('selection_%s.pdf' % x, bbox_inches='tight')

    sfn   = {}

    ## 2D interpolation of HSC selection. 
    f      = interpolate.RegularGridInterpolator((xedges[:-1], yedges[:-1]), HH, bounds_error=False, fill_value=0.0, method='nearest') 
    sfn[x] = f

    return sfn
    

if __name__ == "__main__":
    cm = plt.cm.get_cmap('viridis')

    pd.set_option('display.max_rows',   1000)
    pd.set_option('display.max_columns',  60)
    pd.set_option('display.width',      2000)

    get_SubaruSelectionfn(printit = True, plotit = True)
