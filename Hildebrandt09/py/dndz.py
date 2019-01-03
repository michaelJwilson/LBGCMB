import  numpy  as np
import  pylab  as pl
import  os


def getpz_H09():
  root       = os.environ['LBGCMB']
  data       = np.loadtxt(root + '/Hildebrandt09/dndz/udrops_r23p0t25p5_Masim.dat')

  dz         = data[1,0] - data[0,0]
  norm       = np.sum(data[:,1] * dz)

  pz         = data[:,1] / norm
  midz       = data[:,0] + 0.5 * dz

  return midz, pz


if __name__ == "__main__":
  midz, pz = getpz_H09()

  pl.plot(midz, pz)

  pl.show()
