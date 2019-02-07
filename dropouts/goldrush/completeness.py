import  os
import  numpy              as      np
import  pylab              as      pl

from    scipy.interpolate  import  interp1d 
from    goldrush.specs     import  samplestats
from    elg_nz             import  elg_pz


def interp_completeness(zee, drop='g'):
  stats          =  samplestats(printit=False)

  root           =  os.environ['LBGCMB']
  path           =  root + '/dropouts/goldrush/cats/completeness/completeness_z%d.dat' % round(stats[drop]['z'])

  completeness   =  np.loadtxt(path)
  interp         =  interp1d(completeness[:,0], completeness[:,1], kind='linear', assume_sorted=False, bounds_error=False, fill_value=0.0)

  return  interp(zee)

def plot_completeness():
  zs  = np.arange(2., 8., 0.1)

  for band in ['g', 'r', 'i', 'z']:
    pl.plot(zs, interp_completeness(zs, band), label=band)

  pl.xlabel(r'$z$')
  pl.ylabel(r'Completeness')

  pl.title(r'HSC Goldrush')
  pl.legend(loc=2, ncol=1)  

  pl.show()
  ##  pl.savefig('../plots/completeness.pdf')

def get_dropoutpz(drop='g'):
  ##  Note:  returns normalised completeness curve. 
  ##         requires LF and volume factors to be added. 
  stats          = samplestats(printit = False)

  root           = os.environ['LBGCMB']
  path           = root + '/dropouts/goldrush/cats/completeness/completeness_z%d.dat' % round(stats[drop]['z'])

  data           = np.loadtxt(path)
  
  zee            = data[:,0]
  completeness   = data[:,1]
  
  dz             = zee[1] - zee[0]
  norm           = np.sum(completeness * dz)

  pzee           = completeness / norm

  return  zee, pzee

def calc_moments(zee, pz):
  dz             = zee[1] - zee[0]
  mean_zee       = np.sum(zee * pz * dz)

  var_zee        = np.sum(zee * zee * pz * dz) - mean_zee ** 2.
  std_zee        = np.sqrt(var_zee)

  print(mean_zee, 2. * std_zee, 2. * np.sqrt( 2. * np.log(2.) ) * std_zee)
  

if __name__ == "__main__":
  print("\n\nWelcome to completeness.\n\n")

  
  plot_completeness()

  '''
  for band in ['g', 'r']:
    z, pz = get_dropoutpz(band)

    calc_moments(z, pz)
  '''
  print("\n\nDone.\n\n")
