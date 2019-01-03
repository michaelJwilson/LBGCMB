import  re
import  os
import  time
import  zwarning

import  desi
import  rrio
import  json               
import  collections
import  pylab                as      pl
import  numpy                as      np

from    desispec.resolution  import  Resolution
from    astropy.table        import  Table
from    utils                import  sci_notation


class PlotSpec(object):
    def __init__(self, infile, itarget=0, survey='desi', type='Shapley', Quantile=3, latex=True, printit=True):
        '''
        TO DO: document
        '''

        ##  Isolate imports of optional dependencies
        import  warnings
        import  matplotlib
        import  matplotlib.pyplot  as  plt


        specfile                  = infile
        output                    = infile.split('/')[-1].split('.fits')[0]

        if type == 'Shapley':
          subpath = type + '/Q%d/' % Quantile
        
        else:
          subpath = type  

        rrh5file                  = os.environ['LBGCMB'] + '/redrock/dat/rrh5/%s/%s/%s.h5'          % (survey, subpath, output)
        zbestfile                 = os.environ['LBGCMB'] + '/redrock/dat/zbest/%s/%s/%s.fits'       % (survey, subpath, output)
        metafile                  = os.environ['LBGCMB'] + '/redrock/dat/meta/%s/%s/%s_meta.json'   % (survey, subpath, output)                          
        truth_path                = os.environ['LBGCMB'] + '/quickspectra/dat/in_quickspectra/%s/'  %  type
        
        self.survey               = survey
        self.targets, targetids   = desi.read_spectra([specfile])
    
        self.templates            = rrio.read_templates()
        self.zscan, self.zfit     = rrio.read_zscan(rrh5file)
        self.zbest                = Table.read(zbestfile)  
        self.itarget              = itarget                     
        self.znum                 = 0                               ## Rank of best-fitting redshifts.   
        self.smooth               = 1
        self.rasterized           = True
        self.latex                = latex 
        self.metadata             = json.load(open(metafile), object_pairs_hook=collections.OrderedDict)
        self.fluxtruth            = truth_path + self.metadata['SPECTRAFILES'][0].split('/')[-1][:-5].split('_exp')[0] + '.dat'

        print('\n\nReading:  \n%s\n%s\n%s\n%s\n%s\n\n' % (specfile, rrh5file, zbestfile, metafile, self.fluxtruth))
        
        if survey == 'pfs':
          self.dof                = 9300                            ## Extended wavelength coverage. 

        elif survey == 'desi':
          self.dof                = 6896  

        else:
          raise ValueWarning('Degrees of freedom for chi sq. is not available for %s.' % survey)  

        warnings.warn('DOF hardcoded in plotspec.py;  Should match output from redrock.', UserWarning)

        redshift                  = rrh5file.split('z')[1][:3]
        exposure                  = rrh5file.split('exp')[1][:-3]

        self.title                = ''

        if type == 'BC03':
          ## Magnitudes for BC03 production run;  Targets ordered by magnitude.  
          magnitudes              = 20.0 + np.linspace(0., 6.00, 25)
        
          if self.latex: 
            self.title            = r'BC03 for $z$' + '=%.2lf' % np.float(redshift) + '; Exposure: %ss' % str(exposure) 
            self.title           += '; Magnitude: %.2lf;  '    % magnitudes[self.itarget]

        else:
            if self.latex:
              self.title          = r'Q3 for $z$' + '=%.2lf' % np.float(redshift) + '; Exposure: %ss;  ' % str(exposure)    
        
        self._fig                 = plt.figure()

        ntemplates                = len(self.templates)
        
        self._ax                  = []
        
        ## Chi2 plots on top of plotting spectra against all available best-fitting templates. 
        self._ax.append(self._fig.add_subplot(1 + ntemplates, 1, 1))
        
        for n in np.arange(ntemplates):
            self._ax.append(self._fig.add_subplot(1 + ntemplates, 1, 2 + n))
        
        print('\n\nPlotting redrock output for TARGETID: %s' % self.itarget)
        
        self.chi2_plot()
        self.spectrum_plot()
        
    def chi2_plot(self, keepzoom=False, interactive = False, title='', plot_bestfit=True):
        from    scipy.signal        import medfilt
        from    collections         import OrderedDict

        import  matplotlib
        import  matplotlib.pyplot   as     plt
        

        target     = self.targets[self.itarget]    
        zfit       =    self.zfit[self.zfit['targetid'] == target.id]
        
        self.nznum = len(zfit)
        
        zz         = zfit[zfit['znum'] == self.znum][0]
        coeff      = zz['coeff']
        
        for spectype in self.zscan[target.id]:
            print  'Spectype available: %s' % spectype
        
        for ax in self._ax:
            ax.clear()
        
        colors      = ['k', 'b', 'c', 'g', 'r']

        bestfitChi2 = zfit['chi2'][zfit['chi2'] == zfit['chi2'].min()] 
        
        for i, spectype in enumerate(self.zscan[target.id]):
            label       = re.escape(str(spectype))
            zscaling    = 1.0

            zx          = self.zscan[target.id][spectype]
                
            ## if spectype.split(':')[0].upper() == 'STAR':
            ##   zscaling  = 1.e0
            ##   label    += ' (z scaling: %.0lf)' % zscaling
            
            self._ax[0].plot(zscaling * zx['redshifts'], zx['zchi2'] / self.dof, colors[i] + '-', rasterized=self.rasterized, alpha=1.0, label=label)
        
        if plot_bestfit:
          ymin  = self._ax[0].get_ylim()[0]
          ymax  = self._ax[0].get_ylim()[1]

          ytext = ymin + 0.6*(ymax - ymin)

          for row in zfit: 
           self._ax[0].plot(row['z'], row['chi2'] / self.dof, 'r.', label='', markersize = 4. + row['znum'], rasterized=False, markerfacecolor='None', lw=.1)

          ## for row in zfit:
          ##   self._ax[0].text(row['z'], ytext, str(row['znum']), verticalalignment='top')
        
        err_str = sci_notation(zz['zerr'], decimal_digits=2, precision=None, exponent=-3)

        if self.latex:
          self.title  =  '%s' % self.survey.upper() + ' Target %d:  ' % target.id + self.title

          self.title +=  'Best $z$: %.3lf ' % zz['z']
          self.title += r' $\pm$'  + ' %s ' % err_str
          self.title +=  ' (%s)' % re.escape(zz['spectype'])        
        
        self._ax[0].legend(loc=1, ncol=1)
        self._ax[0].set_title(self.title)
        
        self._ax[0].set_ylabel(r'$\chi_{\nu}^2$')
        self._ax[0].set_xlabel(r'$z$')
                
        self._ax[0].xaxis.labelpad = -2.0
        
    def spectrum_plot(self, keepzoom=False, interactive = False, title=''):   
        from    scipy.signal        import medfilt

        target     = self.targets[self.itarget]
        zfit       =    self.zfit[self.zfit['targetid'] == target.id]

        for n, ax in enumerate(self._ax[1:]):            
            ax.clear()

            ymin   = ymax = 0.0

            tp     = self.templates[n] 
            
            zz     = zfit[zfit['spectype'] == tp.type][0]
            coeff  = zz['coeff']
            
            for spec in target.coadd:
                print 'Number of coadds: %d for target: %d' % (len(target.coadd), self.itarget)
                
                mx    =  tp.eval(coeff[0:tp.nbasis], spec.wave, zz['z']) * (1. + zz['z'])
                model =  spec.R.dot(mx)

                flux  =  spec.flux.copy()

                isbad = (spec.ivar == 0)
                
                mx[isbad]    = np.NaN
                flux[isbad]  = np.NaN
                model[isbad] = np.NaN
        
                ax.plot(spec.wave, np.tanh(medfilt(flux,  self.smooth)),   'k-',    label='Flux', zorder=1, rasterized=self.rasterized)

                ## Smoothed by instrument resolution. 
                ax.plot(spec.wave, np.tanh(medfilt(mx,    self.smooth)), c='gold',  label='Model * R', zorder=3, alpha=0.7, rasterized=self.rasterized)  

                ymin = min(ymin, flux[~isbad].min())
                ymax = max(ymax, flux[~isbad].max() * 1.1)
            
            ## Label object type and redshift
            label    = '{}, $z$ = {:.3f}'.format(re.escape(tp.fulltype), zz['z']) 
            label   += r'$\ \pm \ $' 
            label   += sci_notation(zz['zerr'], decimal_digits=2, precision=None, exponent=-3)
            label   += r' $ (\chi_{\nu}^2 = $' + '%.3lf$ )$' % (np.array(zz['zzchi2']).min() / self.dof)

            print('target {} id {} {}'.format(self.itarget, target.id, label))
            
            ymin  =        ax.get_ylim()[0]
            ymax  =  0.9 * ax.get_ylim()[1]

            ytext = ymin + 0.1*(ymax - ymin)

            ax.text(7500., ytext, label, bbox=dict(facecolor='white', alpha=0.9))
            
            if self.fluxtruth is not None:
                ## Plot as simple spectra. 
                spectra = np.loadtxt(self.fluxtruth)

                ax.plot(spectra[:,0], np.tanh(spectra[:, self.itarget + 1]), c='deepskyblue', zorder=2, rasterized=self.rasterized)

            ax.set_ylim([ymin, ymax])
            '''
            ## ZWARN labels
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
                
            ax.text(4.5e3, ytext, label, horizontalalignment='left', color=color, bbox=dict(facecolor='white', alpha=0.9))
            '''
            ax.axhline(0., color='k', alpha=0.4)

            if self.survey == 'pfs':
                ax.set_xlim([3.5e3, 1.3e4])
                ax.set_ylim([-1.05, 1.050])
            
            elif self.survey == 'desi':
                ax.set_xlim([3.6e3, 1.0e4])
                ax.set_ylim([-1.05, 1.050])

            else:
                raise ValueError('\n\nSpectral coverage undefined for %s survey.' % survey)

            ax.set_ylabel(r'$T(F_{\lambda} / 10^{-17} / \ \rm{ergs}/s/\rm{cm}^2/\AA$)')

            if n == (len(self._ax[1:]) - 1):
              ax.set_xlabel(r'Wavelength [$\AA$]')            

            else:
              ax.xaxis.set_visible(False)
 
            ## handles, labels = ax.get_legend_handles_labels()
            ## by_label        = OrderedDict(zip(labels, handles))

            ## leg = ax.legend(by_label.values(), by_label.keys(), loc=1)

            ## leg.get_frame().set_alpha(1.0)
        
