import  os
import  pylab             as      pl
import  numpy             as      np

from    plotspec          import  PlotSpec
from    astropy.table     import  Table
from    utils             import  latexify


latexify(columns=2, equal=False, fontsize=6, ratio=1.25, ggplot=True, usetex=True)

if __name__ == "__main__":    
  print("\n\nWelcome to plot redrock.\n\n")


  survey     =  'desi'   
  ## survey  =  'pfs'
  
  ## type    =  'BC03'
  ## subpath =  ''

  type       =  'Shapley'
  quantile   =         3                   ##  Shapley composite quantile number. 
  subpath    =   type + '/Q%d/' % quantile

  itarget    =         0                   ##  Ordered by magnitude. 
  
  redshift   =       2.0

  exposure   =  60. * 15. * 1              ##  Number of 15 minutes exposures   
  
  ## infile  = os.environ['LBGCMB'] + '/quickspectra/dat/quickspectra/%s/%s/spec-BC03-z%.1lf_exp%d.fits' % (survey, type, redshift, exposure)
  infile     = os.environ['LBGCMB'] + '/quickspectra/dat/quickspectra/%s/%s/spec-Shapley-Q3_199-z%.1lf_exp%d.fits' % (survey, subpath, redshift, exposure)  

  output     = infile.split('/')[-1].split('.fits')[0]
  
  zbestfile  = os.environ['LBGCMB'] + '/redrock/dat/zbest/%s/%s/%s.fits' % (survey, subpath, output)

  print
  print(Table.read(zbestfile))
  print
  
  p          = PlotSpec(infile, itarget = itarget, survey=survey, type=type, Quantile = quantile, latex=True)
  
  pl.savefig(os.environ['LBGCMB'] + '/redrock/plots/%s/rrplotspec.pdf' % type, bbox_inches='tight')        
  
  print("\n\nDone.\n\n")
