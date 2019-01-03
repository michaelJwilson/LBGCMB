import  os
import  rrio
import  pylab              as      pl
import  numpy              as      np
import  matplotlib.pyplot  as      plt

from    plotspec           import  PlotSpec
from    astropy.table      import  Table
from    utils              import  latexify


if __name__ == "__main__":
  print("\n\nWelcome to exposure.\n\n")

  
  survey          =  'desi'
  sample_rate     =      1

  ##  BC03 production run.                                                                                                          
  ##  redshifts   =  3.5 + np.linspace(0., 2.50, 25)[::sample_rate]
  ##  magnitudes  = 20.0 + np.linspace(0., 6.00, 25)[::sample_rate]

  ##  Shapley production run. 
  redshifts       =   1.5 + np.linspace(0., 1.50, 15)[::sample_rate]
  magnitudes      =  20.0 +   np.arange(0., 6.00,  1)

  exposures       =   60. * 15. * np.arange(1, 30, 3)[::sample_rate]

  print('\n\nSolving for redshifts:\n'  + ''.join('%.3lf, ' % x for x in redshifts))
  print('\n\nSolving for magnitudes:\n' + ''.join('%.3lf, ' % x for x in magnitudes))
  print('\n\nSolving for exposures:\n'  + ''.join('%.3lf, ' % x for x in exposures))

  ##  type        =  'BC03'                                                                                                                                 
  ##  subpath     =  ''                                                                                                                                     

  ##  Shapley quantile number.                                                                                                                               
  type            =  'Shapley'
  Quantile        =         3
  
  ##  subpath     =  type + '/Q%d' % Quantile
  subpath         =  type + '/all/'

  template_types  =  ['Q4_198', 'QX_811', 'Q3_199', 'Q2_198', 'Q0_199']  
  result          =  {key: [] for key in template_types}
  
  for redshift in redshifts:
    for exposure in exposures:
      ##  infile   =  os.environ['LBGCMB'] + '/quickspectra/dat/BC03/%s/LBGs/spec-BC03-z%.1lf_exp%d.fits'         % (survey,   redshift, exposure)
      ##  infile   =  os.environ['LBGCMB'] + '/quickspectra/dat/Shapley/Q3/spec-Shapley-Q3_199-z%.3lf_exp%d.fits' % (redshift, exposure)
      infile       =  os.environ['LBGCMB'] + '/quickspectra/dat/Shapley/all/spec-Shapley-all-z%.3lf_exp%d.fits'   % (redshift, exposure)  

      output       =  infile.split('/')[-1].split('.fits')[0]
        
      zbestfile    =  os.environ['LBGCMB'] + '/redrock/dat/zbest/%s/%s/%s.fits' % (survey, subpath, output)
      rrh5file     =  os.environ['LBGCMB'] + '/redrock/dat/rrh5/%s/%s/%s.h5'    % (survey, subpath, output)

      Data         =  Table.read(zbestfile)

      ## print
      ## print(Data)
      ## print

      zscan, zfit  =  rrio.read_zscan(rrh5file)
      
      for kk, template_type in enumerate(template_types):
        ##  Target ID ordered by magnitude above. 
        for jj, magnitude in enumerate(magnitudes):
          id          = jj + kk * len(magnitudes)
          lbg_rchi2   = zscan[id][template_type]['zchi2'].min()
        
          ##  Rank of best fitting redshift.
          znum        = 0 
          zz          = zfit[zfit['znum'] == 0][0]

          ##  Find if redrock warning on redshift, or if fitted redshift is 5 sigma from truth. 
          if (zz['zwarn'] != 0):
            success   = 0

            print('\nCaught by ZWARN: z= %.3lf, m=%.3lf, exp=%d; z hat = %.3lf +- %.3le, X2 = %.3lf, warning: %s, success: %d' % (redshift, magnitude,\
                                                                                                                                  exposure, zz['z'],\
                                                                                                                                  zz['zerr'], lbg_rchi2,\
                                                                                                                                  zz['zwarn'], success))
          
          elif  (np.abs((redshift - zz['z'])) > 0.02) & (zz['zwarn'] == 0):
          ## elif  (np.abs((redshift - zz['z'])) > 5. * zz['zerr']) & (zz['zwarn'] == 0):  
            success   = 0

            print('\nCaught by ZOUTLIER: z= %.3lf, m=%.3lf, exp=%d; z hat = %.3lf +- %.3le, X2 = %.3lf, warning: %d, success: %d' % (redshift, magnitude,\
                                                                                                                                     exposure, zz['z'],\
                                                                                                                                     zz['zerr'], lbg_rchi2,\
                                                                                                                                     zz['zwarn'], success))
          
          else:
            success   = 1

          ##  Now check on successful discrimination from other templates.   
          for k, spectype in enumerate(zscan[id]):
            ##  Exclude all Shapley composites.  i.e. do not rule out a 'successful' redshift if Q0 is an equally good fit to Q4 fit. 
            if not spectype in template_types:
              zx      = zscan[id][spectype]
              rchi2   = zx['zchi2'].min()

              ##  Difference in chi^2 must be at least twenty.
              if (rchi2 - lbg_rchi2) > 20:
                continue

              else:
                print('\nCaught by TEMPLATE DEGENERACY: z= %.3lf, m=%.3lf, exp=%d; z hat = %.3lf +- %.3le, X2 = %.3lf, template: %s, confusion: %s at z=%.3lf' % (redshift,\
                                                                                                                                                                  magnitude,\
                                                                                                                                                    exposure, zz['z'],\
                                                                                                                                                    zz['zerr'], lbg_rchi2,\
                                                                                                                                                    template_type, spectype,\
                                                                                                                                                    zx['zfit']['z'][0]))
                
                success = 0 

          if success == 1:
            ##  Save only those with successfully identified redshift. 
            result[template_type].append([redshift, magnitude, exposure, success])

  print('\n\n\nExposure times required:\n\n')

  for template_type in template_types:
    print('\n\nFor template:  %s' % template_type)

    rresult    = np.array(result[template_type])
    rrresult   = []

    ##  Loop over the unique redshifts.  Note:  Return of np.unique is sorted. 
    for redshift in np.unique(rresult[:,0]):
      for magnitude in np.unique(rresult[rresult[:,0] == redshift][:, 1]):
        slice   = rresult[rresult[:,0] ==  redshift] 
        sslice  =   slice[slice[:,1]   == magnitude]
      
        ## print
        ## print redshift, magnitude
        ## print sslice

        ##  Append shortest exposure time with successful redshift for given z and mag.
        min_exptime = sslice[:,2].min()

        rrresult.append([redshift, magnitude, min_exptime])

        print('z=%.3lf, m=%.3lf, min. exposure time is %.3lf.' % (redshift, magnitude, min_exptime))

    
    rrresult = np.array(rrresult) 

    ##  Grid of redshifts and magnitudes. 
    xv, yv  = np.meshgrid(redshifts, magnitudes)

    ##  Load with np.nan (coloured white) where no successful redshift is available. 
    zv      = np.zeros_like(yv) * np.nan

    for row in rrresult:
      ##  Where a successful redshift was available, load with shortest exposure time. 
      zv[(xv == row[0]) & (yv == row[1])] = row[2] / 60. / 60.

    zv      = np.ma.masked_invalid(zv)
    
    pl.clf()

    ##  And plot.
    latexify(fig_width=None, fig_height=None, columns=1, equal=True, ggplot=False)

    fig     = pl.gcf()

    cmap    = plt.get_cmap('viridis')
    cmap.set_bad(color = 'white', alpha = 1.)

    plt.pcolormesh(xv, yv, zv, cmap = cmap, vmin = 0., vmax = 6.)

    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Exposure time [hours]', rotation=270, labelpad=20)

    pl.xlim(redshifts[0], redshifts[-1])

    pl.xlabel(r'$z$')
    pl.ylabel(r'$g_{\rm{AB}}$')

    ## pl.savefig(os.environ['LBGCMB'] + '/redrock/plots/%s/exposure_%s.pdf' % (type, survey), bbox_inches='tight')
    pl.savefig(os.environ['LBGCMB'] + '/redrock/plots/tmp/exposure_%s_%s.pdf' % (survey, template_type), bbox_inches='tight')
    
  print("\n\nDone.\n\n")
