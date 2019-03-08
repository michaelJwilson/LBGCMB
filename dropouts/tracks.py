import  json
import  glob
import  matplotlib         as      mpl
import  numpy              as      np
import  pylab              as      pl

import  matplotlib.cm      as      cm
import  matplotlib.pyplot  as      plt

from    collections        import  OrderedDict
from    pylab              import  rcParams
from    utils              import  latexify


latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=12, ggplot=True, usetex=True)

def colourtrack(dband = 'g'):    
                  ## Hildebrandt (2018), https://arxiv.org/pdf/0903.3951.pdf;  Otherwise, GoldRush (Ono, 2018).  
    hsc        = {'u': {'bcol': 'u-g', 'rcol': 'g-r', 'minbcol': 1.5, 'maxrcol': 1.2, 'gradient': 1.5, 'intercept': 0.75, 'hiz': 4.5},\
                  'g': {'bcol': 'g-r', 'rcol': 'r-i', 'minbcol': 1.0, 'maxrcol': 1.0, 'gradient': 1.0, 'intercept': 1.00, 'hiz': 6.0},\
                  'r': {'bcol': 'r-i', 'rcol': 'i-z', 'minbcol': 1.2, 'maxrcol': 0.7, 'gradient': 1.5, 'intercept': 1.00, 'hiz': 7.5},\
                  'z': {'bcol': 'i-z', 'rcol': 'z-y', 'minbcol': 1.5, 'maxrcol': 0.5, 'gradient': 2.0, 'intercept': 1.10, 'hiz': 7.5}}
    
    bcol       =  hsc[dband]['bcol']
    rcol       =  hsc[dband]['rcol']

    minbcol    =  hsc[dband]['minbcol']        ## Detect the break
    maxrcol    =  hsc[dband]['maxrcol']        ## Flat spectra above break. 
    
    gradient   =  hsc[dband]['gradient']
    intercept  =  hsc[dband]['intercept']

    rcols      =  np.arange(-1.0, 4.0, 0.001)
    bcols      =  gradient * rcols + intercept

    minrcol    = (minbcol - intercept) / gradient
    bcollim    =  gradient * maxrcol + intercept

    ## rcParams['figure.figsize'] = 3.5, 3.5
    
    pl.clf()

    ## Plot selection box perimeter. 
    pl.plot(rcols[rcols < minrcol], minbcol * np.ones_like(rcols[rcols < minrcol]), 'k', lw=0.4)
    pl.plot(maxrcol * np.ones_like(bcols[bcols > bcollim]), bcols[bcols > bcollim], 'k', lw=0.4)
    
    ## Gradient
    pl.plot(rcols[(rcols > minrcol) & (bcols < bcollim)], bcols[(rcols > minrcol) & (bcols < bcollim)], c='k', linestyle='-', lw=0.4)

    ## Stars.                                                                                                                                             
    stars   = []

    for file in glob.glob("./stellar_spectra/colors/*.json"):
      colors   = json.load(open(file), object_pairs_hook = OrderedDict)
      stars.append([0.0, colors[rcol], colors[bcol]])

    stars  = np.array(stars)

    ## print(stars)                                                                                                                                                                                                     
    plt.scatter(stars[:,1], stars[:,2], c='k', s=2)

    ## LRG tracks; 
    reds   = plt.get_cmap('Reds')
    tracks = []

    for file in glob.glob("./colors/LRG/*.json"):
        colors   = json.load(open(file), object_pairs_hook = OrderedDict)

        string   = file.split('z')[1]
        string   = string[:4]

        redshift = np.float(string)

        tracks.append([redshift, colors[rcol], colors[bcol]])

    tracks  = np.array(tracks)
    tracks  = tracks[tracks[:,0].argsort()]

    np.set_printoptions(precision=3)

    ## print(tracks)                                                                                                                                         

    tracks  = tracks[tracks[:,0] < hsc[dband]['hiz']]
    cax     = plt.scatter(tracks[:,1], tracks[:,2], c = reds(tracks[:,0] / tracks[:,0].max()), s=9, marker='s')

    ## LBG tracks;  Color ordering: ["u-g", "g-r", "r-i", "i-z"]
    tracks = []

    for file in glob.glob("./colors/LBG/*.json"):
        colors   = json.load(open(file), object_pairs_hook = OrderedDict)

        string   = file.split('z')[1]
        string   = string[:4]

        redshift = np.float(string)
                
        tracks.append([redshift, colors[rcol], colors[bcol]])

    tracks  = np.array(tracks)
    tracks  = tracks[tracks[:,0].argsort()]

    np.set_printoptions(precision=3)

    ## print(tracks)

    tracks  = tracks[tracks[:,0] < hsc[dband]['hiz']]
    cax     = plt.scatter(tracks[:,1], tracks[:,2], c=tracks[:,0], s=9)

    plt.colorbar(cax, label=r'$z$', boundaries = np.linspace(0., 5., 24), format='%.2lf')


    ax = pl.gca()

    ax.set_axis_on()

    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')

    pl.xlabel(r'$' + str(rcol) + '$')
    pl.ylabel(r'$' + str(bcol) + '$')

    pl.xlim(-0.05, 1.2)
    pl.ylim(-0.05, 2.5)

    ## pl.xlim(-1.0, 4.0)
    ## pl.ylim(-1.0, 4.0)

    pl.savefig('plots/tracks/%scolour_track.pdf' % dband, bbox_inches='tight')


if __name__ == "__main__":
  print('\n\nWelcome to tracks.\n\n')
    
  colourtrack(dband = 'u')

  print('\n\nDone.\n\n')
