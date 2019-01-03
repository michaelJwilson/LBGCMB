import pylab as pl
import time
import numpy as np

from   .                   import io
from   .                   import zwarning
from   desispec.resolution import Resolution

class PlotSpec(object):
    def __init__(self, targets, templates, zscan, zfit, truth=None, interactive=False, fluxtruth=None, title=""):
        '''
        TODO: document
        '''

        #- Isolate imports of optional dependencies
        import matplotlib.pyplot as plt

        self.targets   = targets
        self.templates = templates
        self.zscan     = zscan
        self.zfit      = zfit
        self.itarget   = 0              ## Target number
        self.znum      = 0              ## Rank of best-fitting redshifts.   
        self.smooth    = 1
        self.truth     = truth
        self.fluxtruth = fluxtruth

        self._fig      = plt.figure()

        self._fig.set_size_inches(18.5, 10.5)

        self._ax1      = self._fig.add_subplot(311)
        self._ax2      = self._fig.add_subplot(312)
        self._ax3      = self._fig.add_subplot(313)

        if interactive == True:
            self._cid      = self._fig.canvas.mpl_connect('key_press_event', self._onkeypress)

            #- Instructions
            print("-------------------------------------------------------------------------")
            print("Select window then use keyboard shortcuts to navigate:")
            print("    up/down arrow: previous/next target")
            print("    left/right arrow: previous/next redshift fit for this target")
            print("    (d)etails")
            print("-------------------------------------------------------------------------")

            #- Disable some default matplotlib key bindings so that we can use keys
            #- TODO: cache and reset when done
            plt.rcParams['keymap.forward'] = ''
            plt.rcParams['keymap.back']    = ''

            plt.ion()

        self.plot(interactive=interactive, title=title)

    def _onkeypress(self, event):
        if event.key == 'right':
            self.znum = (self.znum + 1) % self.nznum
            self.plot(keepzoom=True)

        elif event.key == 'left':
            self.znum = (self.znum - 1) % self.nznum
            self.plot(keepzoom=True)

        elif event.key == 'down':
            if self.itarget == len(self.targets)-1:
                print('At last target')
            else:
                self.znum = 0
                self.itarget += 1
                self.plot()

        elif event.key == 'up':
            if self.itarget == 0:
                print('Already at first target')
            else:
                self.znum = 0
                self.itarget -= 1
                self.plot()

        elif event.key == 'd':
            target = self.targets[self.itarget]    
            zfit   = self.zfit[self.zfit['targetid'] == target.id]
            
            print('target {}'.format(target.id))
            print(zfit['znum', 'spectype', 'z', 'zerr', 'zwarn', 'chi2'])

    def plot(self, keepzoom=False, interactive = False, title=""):
        from    scipy.signal        import medfilt
        from    collections         import OrderedDict

        import  matplotlib.pyplot   as     plt
        
        target = self.targets[self.itarget]    
        zfit   = self.zfit[self.zfit['targetid'] == target.id]
        
        self.nznum = len(zfit)
        
        zz     = zfit[zfit['znum'] == self.znum][0]
        coeff  = zz['coeff']

        if keepzoom:  ## zscan plot
            force_xlim = self._ax1.get_xlim()
            force_ylim = self._ax1.get_ylim()

        self._ax1.clear()
        
        for spectype, fmt in [('STAR', 'k-'), ('GALAXY', 'b-'), ('QSO', 'g-'), (u'LBG:NOLINES', 'c-')]:
            if spectype in self.zscan[target.id]:
                print "Spectype available: %s" % spectype
                
                zx = self.zscan[target.id][spectype]

                self._ax1.plot(zx['redshifts'], zx['zchi2'], fmt, alpha=0.2, label='_none_')
                self._ax1.plot(zx['redshifts'], zx['zchi2'] + zx['penalty'], fmt, label=spectype)

        self._ax1.plot(zfit['z'], zfit['chi2'], 'r.', label='_none_')

        for row in zfit:
            self._ax1.text(row['z'], row['chi2'], str(row['znum']), verticalalignment='top')

        if self.truth is not None:
            i = np.where(self.truth['targetid'] == target.id)[0]

            if len(i) > 0:
                ztrue = self.truth['ztrue'][i[0]]

                self._ax1.axvline(ztrue, color='g', alpha=0.5)

            else:
                print('WARNING: target id {} not in truth table'.format(target.id))

        self._ax1.axvline(zz['z'],    color='k', alpha=0.4)
        self._ax1.axhline(zz['chi2'], color='k', alpha=0.4)

        self._ax1.legend(loc=2)

        self._ax1.set_title('%s;  ' % title + 'Target {};  Best fitted $z$: {:.3f};  Type: {}'.format(target.id, zz['z'], zz['spectype']))
        self._ax1.set_ylabel(r'$\chi^2$')
        self._ax1.set_xlabel('Redshift')

        self._ax1.set_xlim(0.0, 2.9)

        if keepzoom:
            self._ax1.set_xlim(*force_xlim)
            self._ax1.set_ylim(*force_ylim)
    
        #----- Spectrum plot
        for ax in [self._ax2, self._ax3]:
            if keepzoom:
                force_xlim = ax.get_xlim()
                force_ylim = ax.get_ylim()
            
            ax.clear()

            ymin = ymax = 0.0

            for tp in self.templates:
                if ax == self._ax2 and tp.type == 'LBG':
                    ## zz derived from zfit; info. on best fitting template irrespective of type.  
                    print("\n\nPlotting Galaxy fit.")

                    break

                if ax == self._ax3 and tp.type == 'GALAXY':
                    print("\n\nPlotting LBG fit.")
                    
                    break
                    
            zz     = zfit[zfit['spectype'] == tp.type][0]
            coeff  = zz['coeff']
                
            for spec in target.coadd:
                print "Number of coadds: %d for target: %d" % (len(target.coadd), self.itarget)
                
                mx    =  tp.eval(coeff[0:tp.nbasis], spec.wave, zz['z']) * (1. + zz['z'])
                model =  spec.R.dot(mx)

                flux  =  spec.flux.copy()

                isbad = (spec.ivar == 0)
                
                mx[isbad]    = np.NaN
                flux[isbad]  = np.NaN
                model[isbad] = np.NaN
        
                ax.plot(spec.wave, medfilt(flux,  self.smooth), 'k-',  label='Flux')
                ax.plot(spec.wave, medfilt(mx,    self.smooth), 'r-',  label='Model * R')  ## Smoothed by instrument resolution.  

                ## ymin = min(ymin, np.percentile(flux[~isbad], 1))
                ## ymax = max(ymax, np.percentile(flux[~isbad], 99), np.max(model)*1.05)

                ymin = min(ymin, flux[~isbad].min())
                ymax = max(ymax, flux[~isbad].max() * 1.1)

                ## ymin = min(ymin, model[~isbad].min())                                                                            
                ## ymax = max(ymax, model[~isbad].max() * 1.1)

            ## if self.fluxtruth is not None:
            ##    ax.plot(self.fluxtruth[:,0], self.fluxtruth[:,1])

            #- Label object type and redshift
            label = '{} z={:.3f}'.format(tp.fulltype, zz['z'])
        
            print('target {} id {} {}'.format(self.itarget, target.id, label))
        
            ytext = ymin + 0.9*(ymax - ymin)

            ax.text(6000., ytext, label)

            #- ZWARN labels
            if zz['zwarn'] != 0:
                label = list()
                    
                for name, mask in zwarning.ZWarningMask.flags():
                    if (zz['zwarn'] & mask) != 0:
                        label.append(name)
                            
                label = '\n'.join(label)
                color = 'r'

            else:
                label = 'No warnings'
                color = 'g'

            ax.text(9.85e3, ytext, label, horizontalalignment='right', color=color)

            ax.axhline(0, color='k', alpha=0.2)

            if keepzoom:
                ax.set_xlim(*force_xlim)
                ax.set_ylim(*force_ylim)
        
            else:
                ax.set_ylim(ymin, ymax)
                ax.set_xlim([3.5e3, 1e4])

            ax.set_ylabel('Flux')
            ax.set_xlabel('Wavelength [A]')

            handles, labels = ax.get_legend_handles_labels()
            by_label        = OrderedDict(zip(labels, handles))

            leg = ax.legend(by_label.values(), by_label.keys(), loc=3)

            leg.get_frame().set_alpha(1.0)

            if interactive is True:
                self._fig.canvas.draw()
        
                plt.show()
