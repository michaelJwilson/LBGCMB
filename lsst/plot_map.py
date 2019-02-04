import  json
import  sfdmap
import  sqlite3
import  numpy               as  np
import  pylab               as  pl
import  matplotlib.pyplot   as  plt

from    astropy             import units as u
from    astropy.coordinates import SkyCoord


def radec2project(ra, dec):
    return (np.radians(ra) - np.pi, np.radians(dec))


if __name__ == '__main__':
  with open('lsst_u_fiveyr.json') as json_file:  
    hitmap = json.load(json_file)
                                                                                                                                       
  Area  = 0.0
  YEAR  =   5 

  ##  And plot.                                                                                                                                         
  fig   =  plt.figure(figsize=(8, 8))
  ax    =  plt.subplot(111, projection="aitoff")

  ##  Current expected performance                                                                                                                      
  single_m5 = {'u': 23.98, 'g': 24.91, 'r': 24.42, 'i': 23.97, 'z': 23.38, 'y': 22.47}

  for field in hitmap.keys():
    x, y  = radec2project(hitmap[field]['RA'], hitmap[field]['DEC'])
    cax   = ax.scatter(x, y, c=np.array(hitmap[field]['CoAdd5sig']), cmap=plt.cm.coolwarm, rasterized=True,\
                       vmin = single_m5['u'] - 0.1, vmax = single_m5['u'] + 2.5)

    Area += hitmap[field]['FOV']

  print('\n\nLSST area after Year %d:  %.3lf' % (YEAR, Area))

  plt.grid(True)
  pl.tight_layout()

  pl.colorbar(cax, shrink=0.4)
  
  pl.show()
  ##  pl.savefig('wfd.pdf')
