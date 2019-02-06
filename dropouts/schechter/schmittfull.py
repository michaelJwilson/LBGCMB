import  os
import  numpy              as  np
import  pylab              as  pl
import  matplotlib.pyplot  as  plt


root     = os.environ['LBGCMB']
data     = np.loadtxt(root + '/dat/schmittfull_nbar.dat') 

indices  = np.where(data[:,0] <= 7.0)[0]

data     = data[indices, :]
data     = np.flipud(data)

summed   = np.cumsum(data[:, 1])

pl.plot(data[:,0], summed, 'k-')

ax       = plt.gca()

## ax.invert_xaxis()

pl.savefig("../plots/schmittfull_nbar.pdf")





