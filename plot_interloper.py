import os 
import json
import numpy             as     np
import pylab             as     pl
import matplotlib.pyplot as     plt

from   utils             import latexify
from   fisher_contour    import plot_ellipse


##  And plot contour ...                                                                                                                                                                                                 
pl.clf()

latexify(columns=1, equal=True, fontsize=12)

fig  = plt.gcf()
ax   = plt.gca()

band     =   'u'
int_frac = 0.06 

with open('dat/result4interloper_%s_intfrac_%.2lf.json' % (band, int_frac), 'r') as rfile:
    ##  {'peakz': peakz, 'fid_sig8': np.float(fid_sig8), 'fid_b1': np.float(fid_b1), 'biased_sig8': biased_sig8, 'biased_b1': biased_b1, 'iFisher': iFisher.tolist()}
    data = json.load(rfile)

##  Unpack.
peakz        = data['peakz']
fid_sig8     = data['fid_sig8']
fid_b1       = data['fid_b1']
biased_sig8  = data['biased_sig8']
biased_b1    = data['biased_b1']
iFisher      = np.array(data['iFisher'])

pl.axvline(x=fid_sig8, ymin=0.0, ymax=1.0, c='k', label=r'$(\sigma_8)$', alpha=1., lw=0.4)

for mass_level, color, alpha in zip([0.99, 0.95, 0.68], ['b', 'b', 'b'], [0.2, 0.4, 0.6]):
    plot_ellipse(x_cent = biased_sig8, y_cent = biased_b1, ax = ax, cov = iFisher, mass_level = mass_level,\
                 fill=True, fill_kwargs={'alpha': alpha, 'c': color}, plot_kwargs={'c': color, 'alpha': 0.0})

pl.plot(biased_sig8,    biased_b1,    'w*', markersize=5, label=r'$(\sigma_8, b_1)$')
##  pl.plot(biased_sig8, biased_b1, 'k*', markersize=5, alpha=0.4, label=r'$(\hat \sigma_8, \hat b_1)$')   

##  pl.legend(frameon=False)

pl.xlabel(r'$\sigma_8(z=%.1lf)$' % peakz)
pl.ylabel(r'$b_1(z=%.1lf)$'      % peakz)

plt.tight_layout()

##  pl.show()
pl.savefig('plots/interloper_bias_contours_%s-drops.pdf' % band)  
