## import  matplotlib;                matplotlib.use('Agg')

import  copy
import  json
import  glob
import  numpy              as      np
import  pylab              as      pl

import  matplotlib.cm      as      cm
import  matplotlib.pyplot  as      plt
import  pandas             as      pd

from    collections        import  OrderedDict
from    pylab              import  rcParams
from    app_mags           import  get_colors
from    utils              import  latexify
from    colourcut_dNdz     import  colourcut
from    astropy.table      import  Table


##  latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=False)

def colourtrack(ttype = 'g'):    
    '''
    Plot locus of colour box for given dropout selection. 
    '''
                  ##  Hildebrandt (2018), https://arxiv.org/pdf/0903.3951.pdf;  Otherwise, GoldRush (Ono, 2018).  
                  ##  u-drops:  1.5 < (u-g) and -1 < (g-r) < 1.2 and 1.5 (g-r) < (u-g) - 0.75.
                  ##      BzK:  

    selections = {'u':   {'bcol': 'u-g', 'rcol': 'g-r', 'minbcol': 1.5, 'maxrcol': 1.2, 'gradient': 1.5, 'intercept': 0.75, 'hiz': 4.5},\
                  'g':   {'bcol': 'g-r', 'rcol': 'r-i', 'minbcol': 1.0, 'maxrcol': 1.0, 'gradient': 1.0, 'intercept': 1.00, 'hiz': 6.0},\
                  'r':   {'bcol': 'r-i', 'rcol': 'i-z', 'minbcol': 1.2, 'maxrcol': 0.7, 'gradient': 1.5, 'intercept': 1.00, 'hiz': 7.5},\
                  'z':   {'bcol': 'i-z', 'rcol': 'z-y', 'minbcol': 1.5, 'maxrcol': 0.5, 'gradient': 2.0, 'intercept': 1.10, 'hiz': 7.5}}
            
    if ttype == 'BzK':
        bcols      =  np.arange(-4.0, 7.0, 0.001)
        rcols      =  bcols -0.2

        pl.plot(rcols, bcols, 'k-')
        pl.plot(2.5 * np.ones_like(rcols[rcols > 2.5]), bcols[rcols > 2.5], 'k--')

    elif ttype == 'Euclid':
        pass

    else:
        ##  Dropouts.
        bcol       =  selections[ttype]['bcol']
        rcol       =  selections[ttype]['rcol']
        
        minbcol    =  selections[ttype]['minbcol']        ## Detect the break
        maxrcol    =  selections[ttype]['maxrcol']        ## Flat spectra above break. 
        
        gradient   =  selections[ttype]['gradient']
        intercept  =  selections[ttype]['intercept']

        rcols      =  np.arange(-5.0, 10.0, 0.001)
        bcols      =  gradient * rcols + intercept

        minrcol    = (minbcol - intercept) / gradient
        bcollim    =  gradient * maxrcol + intercept

        pl.plot(rcols[rcols < minrcol], minbcol * np.ones_like(rcols[rcols < minrcol]), 'k', lw=0.6)
        pl.plot(maxrcol * np.ones_like(bcols[bcols > bcollim]), bcols[bcols > bcollim], 'k', lw=0.6)
    
        ##  Gradient
        pl.plot(rcols[(rcols > minrcol) & (bcols < bcollim)], bcols[(rcols > minrcol) & (bcols < bcollim)], c='k', linestyle='-', lw=0.6)

        if ttype == 'u':
            ##  Hildebrandt u-drops has a lower limit on (g-r)
            pl.plot(-1. * np.ones_like(bcols[bcols > minbcol]), bcols[bcols > minbcol], c='k', linestyle='--', lw=0.6)
            

def plot_cat(ttype = 'g', DODEPTH='FULL'):    
    if   ttype == 'u':
         bcol = 'u-g'
         rcol = 'g-r' 

         bmin = -0.5 - 2.5
         bmax =  2.5 + 2.5

         rmin = -0.3 - 2.5
         rmax =  1.2 + 2.5
         
    elif ttype == 'g':
         bcol = 'g-r'
         rcol = 'r-i'
        
         bmin = -0.5 - 2.5
         bmax =  2.5 + 2.5

         rmin = -0.3 - 2.5
         rmax =  1.2 + 2.5

    elif ttype == 'BMBX':
         bcol = 'g-r'
         rcol = 'r-i'

         bmin = -0.5 - 2.5
         bmax =  2.5 + 2.5

         rmin = -0.3 - 2.5
         rmax =  1.2 + 2.5

    elif ttype == 'BzK':
         bcol = 'B-z'
         rcol = 'z-K'

         bmin = -1.0 - 2.5
         bmax =  4.0 + 2.5

         rmin = -0.5 - 2.5 
         rmax =  6.0 + 2.5

    elif ttype == 'Euclid':
         bcol = 'Y-J'
         rcol = 'J-H'

         bmin = -2.0
         bmax =  3.0

         rmin = -2.0
         rmax =  3.0

    else:
        raise ValueError("\n\nSelection type %s is not available with cat_tracks.\n\n" % ttype)

    ##  -- Start -- 
    files         =  glob.glob('colors/scratch/*.fits')
    files         = [x.split('.fits')[0] for x in files]

    files         = [x.split('colors/scratch/3DHST_')[1] for x in files]

    fields        = [x.split('_')[0] for x in files]
    ufields       = list(set(fields))
    nfields       = len(ufields)

    depths        = [xx.split('_')[1].split('_')[0] for xx in files]
    depths        = ['FULL' if depth != 'DefaultDepths' else depth.upper() for depth in depths]
    udepths       =  set(depths)

    count         =  0

    ##  print(fields, depths)

    for ii, xx in enumerate(files):
        depth  = depths[ii]
        field  = fields[ii]

        drops  = []

        if (field == 'UVUDF') & (depth == DODEPTH):
            fname = 'colors/scratch/3DHST_' + xx + '.fits'

            print('Solving for %s' % fname)

            magt  = Table.read(fname, format='fits')
            cols  = magt.colnames

            mags  = copy.copy(cols)
            mags.remove('id')
            mags.remove('zpeak')

            count += 1

            for row in magt:
              if (row['zpeak'] > -99.):
                rowmags = [row[x] for x in mags]
                magdict = dict(zip(mags, rowmags))

                colors  = get_colors(magdict, get_colors=['g-r', 'r-i', 'u-g', 'g-r', 'z-K', 'u-z', 'B-z', 'J-H', 'Y-J'], fname = None)
                is_lbg  =  colourcut(magdict,  dropband=ttype, good=True, fourthlimit=False, BzK_type='all')

                '''
                if magdict['u'] > 98.:
                    ##  Mapped zero flux in template to mag = 99.
                    ##  Do something sensible with the colours.
                    colors['u-g'] = 2.49
                    
                if magdict['g'] > 98.:
                    colors['g-r'] = 2.49
                '''

                if ttype == 'BzK':
                  BzK     = colors['z-K'] - colors['B-z']

                  print('%+.4lf \t %+.4lf \t %+.4lf \t %+.4lf \t %s' % (magdict['B'], magdict['z'], magdict['K'], BzK, str(is_lbg)))

                else:
                  print('%+.4lf \t %+.4lf \t %+.3lf \t %+.3lf \t %s' % (magdict['u'], magdict['g'], magdict['r'], row['zpeak'], str(is_lbg)))
                
                if is_lbg:
                    cax = plt.scatter(colors[rcol], colors[bcol], c=row['zpeak'], s=20, vmin=0.0, vmax=5.0, rasterized=True, alpha=0.8)

                else:
                    cax = plt.scatter(colors[rcol], colors[bcol], c=row['zpeak'], marker='x', alpha=0.3, s=20, vmin=0.0, vmax=5.0, rasterized=True)

            if count > 3:
                break

    ax     = pl.gca()
    fig    = pl.gcf()

    ##  [left, bottom, width, height].
    ##  cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 
    ##  cb     = plt.colorbar(ax=ax, cax = cbaxes, label=r'redshift')  

    if ttype in ['u', 'g', 'BzK', 'Euclid']:
      plt.colorbar(cax, label=r'redshift')

    return  bcol, rcol, [[rmin, rmax], [bmin, bmax]]


if __name__ == "__main__":
  print("\n\nWelcome to colour tracks.\n\n")

  for DODEPTH in ['FULL', 'DEFAULTDEPTHS']:    ##  ['FULL', 'DEFAULTDEPTHS']
    for ttype in ['u', 'g', 'BzK', 'Euclid']:  ##  ['u', 'g', 'BzK', 'Euclid']
          pl.clf()

          bcol, rcol, [[rmin, rmax], [bmin, bmax]] = plot_cat(ttype=ttype, DODEPTH = DODEPTH)  

          colourtrack(ttype=ttype)
          
          pl.xlabel(r'$%s$' % rcol)
          pl.ylabel(r'$%s$' % bcol)

          pl.xlim(rmin, rmax)
          pl.ylim(bmin, bmax)

          ##  plt.tight_layout()

          ##  pl.show(block=True)
          pl.savefig('plots/ccplots/%s_colour_track_%s.pdf' % (ttype, DODEPTH)) ## bbox_inches='tight'

  print("\n\nDone.\n\n")
