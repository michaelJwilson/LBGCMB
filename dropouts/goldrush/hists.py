import numpy             as np
import pylab             as pl
import pandas            as pd
import matplotlib.pyplot as plt


plt.style.use('ggplot')

## (1) Object ID.
## (2) Right ascension.
## (3) Declination.
## (4) Spectroscopic redshift.
## (5) Apparent magnitude.
## (6) UV absolute magnitude.
## (7) The dropout sample in which the source is selected:
##         1 = g-dropout,
##         2 = r-dropout,
##         3 = i-dropout,
##         4 = z-dropout.
## (8) Galaxy/AGN flag (1 = galaxy; 2 = AGN).
## (9) Reference of spectroscopic redshift.

data              = pd.read_csv('../dat/grush_specs.txt', names=['ID', 'ra', 'dec', 'z', 'm', 'MUV', 'sample', 'GalQSO', 'Ref'], header=37, sep='&')

## Cut on only:  H17 = R. Higuchi et al. in preparation.
data              = data[data.Ref == ' This Study  ']

bins              = 10
## bins           = [np.arange(3.0, 6.25, 0.25), np.arange(22.5, 26.75, 0.25)]             
H, xedges, yedges = np.histogram2d(data['z'].values, data['m'].values, bins=bins)

fig               = plt.figure(figsize=(6,6))
ax                = fig.add_subplot(111, title='', aspect='equal')
X, Y              = np.meshgrid(xedges[:-1], yedges[:-1])

plt.pcolormesh(X, Y, H)

plt.colorbar()

pl.xlabel(r'$z$')
pl.ylabel(r'$m$')

pl.savefig('../plots/mz.pdf', bbox_inches='tight')
