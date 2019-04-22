import  matplotlib         as      mpl 
import  numpy              as      np
import  matplotlib.pyplot  as      plt
import  pylab              as      pl

from    utils              import  latexify
from    utils              import  comoving_distance
from    params             import  get_params 


latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=12)

params   = get_params()

############################
ilist    = [10, 11, 12, 13, 14]

drops    = {3.0: [12.33, 4.0, 0.40], 4.0: [12.83, 6.15, 0.45]}

zee      = 3.0
chi_star = comoving_distance(zee)

a        = 1. / (1. + zee) 

mass     = drops[zee][0]
bias     = drops[zee][1]
unknown  = drops[zee][2]

dm       = np.loadtxt("../dat/dm/dm_%.4lf_10.pkr" % a)

for irun in ilist[1:]:
  dm += np.loadtxt(("../dat/dm/dm_%.4lf_{:d}.pkr" % a).format(irun))

dm      /= float(len(ilist))
dm[:,1] *= 2*np.pi**2/dm[:,0]**3

pkr = np.loadtxt("../dat/summary/drop_%.4lf_%.2lf_%.2lf.pkr" % (a, mass, unknown))
pxr = np.loadtxt("../dat/summary/drop_%.4lf_%.2lf_%.2lf.pxr" % (a, mass, unknown))

if  zee  == 3.0:
   labels = ['$P_{gg}$', '$bP_{gm}$', '$b^2P_{mm}$',          '',         '',                '']

elif zee == 4.0:
   labels = [        '',          '',            '',  '$P_{gm}$', '$P_{mm}$', r'$P_{\rm{lin}}$'] 

else:
  raise  ValueError('\n\nSomething went wrong with zee and labels.\n\n')

pl.loglog(pkr[:,0],        pkr[:,1], label=labels[0])
pl.loglog(pxr[:,0],   pxr[:,1]*bias, label=labels[1])
pl.loglog(dm[:,0],  dm[:,1]*bias**2, label=labels[2])

pl.plot(pxr[:,0], pxr[:,1], label=labels[3], dashes=[3,1])
pl.plot(dm[:,0],  dm[:,1],  label=labels[4], dashes=[3,1])

##  Plot linear. 
linear  = np.loadtxt('../dat/thy/pklin_%.4lf.txt' % a)
pl.plot(linear[:,0], linear[:,1], label=labels[5])

pl.legend(loc=1, ncol=1)

pl.xlabel(r'$k\quad [h\,{\rm Mpc}^{-1}]$')
pl.ylabel(r'$P(k)\quad [h^{-3}\,{\rm Mpc}^3]$')

pl.xlim(0.04, 1.)
pl.ylim(1.e1, 1.e5)

##  Add upper axis.                                                                                                                                                                            
ax  = pl.gca()

ax.set_axis_on()

ax.spines['bottom'].set_color('black')
ax.spines['top'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['right'].set_color('black')

ax2 = ax.twiny()
ax2.grid(False)

##  Set new lims. 
ax2.set_xlim(ax.get_xlim())
ax2.set_xscale('log')

new_tickloc = np.array([.04, .1, .3, .6])
new_ticks   = ['%.0lf' % x for x in new_tickloc * chi_star]

ax2.set_xticks(new_tickloc)                                                                                                                                                         
ax2.set_xticklabels(new_ticks)
                                                                                                                                  
ax2.set_xlabel(r"$L = k \cdot \chi_*$")

ax2.set_axis_on()

ax2.spines['bottom'].set_color('black')
ax2.spines['top'].set_color('black')
ax2.spines['left'].set_color('black')
ax2.spines['right'].set_color('black')

##  Finally, add redshift label.
text = plt.text(0.045, 20, r'$z$ = ' + '%.1lf' % zee)
text.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='white'))

pl.savefig('../plots/pk_z%s.pdf' % ('%.2lf' % zee).replace('.', 'p'), bbox_inches='tight')
