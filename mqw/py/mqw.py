import  numpy              as      np
import  matplotlib.pyplot  as      plt

from    prep_Llls          import  prep_Llls
from    pmh                import  invG_bz, linz_bz, Pmm, get_PkInterps
from    prep_camb          import  CAMB
from    utils              import  comoving_distance
from    collections        import  OrderedDict
from    pz2nbar            import  nbar_convert
from    whitebook_pz       import  whitebook_pz
from    utils              import  latexify, sci_notation
from    scipy.interpolate  import  interp1d
from    whitebook_pz       import  const_pz
from    Cij                import  Cij
from    Cgg                import  Cgg

        
def Fisher(Pk_interps, Llls, tNs, tNp, ps, _bs, pp, _bp, dz = 0.1, fsky=0.1, zmin=3.0, zmax=4.0, fover=0.0, percentiles = [], intlp_zs = [], printit=False):
    ##  Returns the variance on \hat Np = F^{-1}_ii.    
    zs      =  np.arange(zmin + dz / 2., zmax + dz / 2., dz)
    result  =  OrderedDict()
    
    for i, zee in enumerate(zs):
      ##  Returns linear bias of spectroscopic galaxies at each redshift. 
      bs    = _bs(zee) 

      ##  Returns linear bias of photometric galaxies at each redshift.                                                                                    
      bp    = _bp(zee)

      ##  Returns number of spectroscopic redshifts in each shell given nbar in per sq. deg.
      Ns    =  nbar_convert(tNs, unit='str') * ps(zee) * dz

      ##  Returns number of photometric redshifts in each shell.
      Np   =   nbar_convert(tNp, unit='str') * pp(zee) * dz
      
      ##  Returns shot noise: wp = Np for fsat = 0.0;  Eqn. (11) of MQW16. 
      wp   =   Np
      ws   =   Ns

      wps  =   fover * np.min([wp, ws])

      '''
      cgg  =   Cij(Pk_interps, Llls, zee - dz / 2.,  zee + dz / 2.)
      '''

      ub   =   lambda z:  1.0
      cgg  =   Cgg(Pk_interps, Llls, zmin, zmax, pp, ub, survey_pz2=ps, bz2=ub, zeff=False)

      ##  Eqn. (44); sum over k has implicit sum over ell and m. 
      result[zee] = {'Ns': Ns, 'Np': Np, 'bp': bp, 'bs': bs, 'cgg': cgg, 'wp': wp, 'num': (bp * bs * cgg) ** 2., 'denom': cgg * (bp * Np)**2. + wp}

      ##  Eqns. (22) - (25) of McQuinn and White. 
      result[zee]['A0i']   =  bp * Np * bs * Ns * cgg + wps  ## < p   s_i >_L
      result[zee]['Aii']   =    (bs * Ns) ** 2. * cgg +  ws  ## < s_i s_j >_L  
      result[zee]['A0i_j'] =       bp * bs * Ns * cgg        ## < p   s_i >_L, j 
      
      ##  Eqn. (43), to be later normalised. 
      result[zee]['beta']  = (Np * bp) ** 2. * cgg 

    ##  Total spec. zs over the whole redshift range. 
    nspec   = 0.0

    ##  ... and normalisation over redshift for beta.                                                                                                    
    bnorm   = 0.0

    A00     = 0.0   ## < p * p >_L
    sshift  = 0.0

    for zee in result:
        nspec  += result[zee]['Ns']

        ##  Array for each L value.
        A00    += result[zee]['denom'] 
        bnorm  += result[zee]['beta']
        sshift += result[zee]['A0i'] ** 2. / result[zee]['Aii'] 
    
    for zee in result: 
        ##  Schur-Limber limit
        result[zee]['beta']         /=  bnorm

        ##  Fractional error from eqn. (44);  Schur-Limber limit of small r (due to shotnoise, or redshift overlap). 
        result[zee]['ratio']         =  result[zee]['num'] / A00

        ##  Diagonal by definition in this limit. 
        result[zee]['Fii']           =  result[zee]['Ns'] * fsky * np.sum((2. * Llls + 1) * result[zee]['ratio'])
        result[zee]['var_ii']        =    1. / result[zee]['Fii']
        result[zee]['ferr_ii']       = np.sqrt(result[zee]['var_ii']) / result[zee]['Np']
    
    ##  Limber approximation, but outwith the Schur limit.    
    for zi in result:
        ##  Schur(L)
        result[zi]['S']              = A00 / (A00 - sshift)
        result[zi]['r']              = result[zi]['A0i'] / np.sqrt(A00 * result[zi]['Aii'])
        
    ##  Non-diagonal Fisher matrix outwith the Schur limit. 
    Fisher = np.zeros((len(zs), len(zs)))

    for i, zi in enumerate(result):
      Fisher[i, i] = fsky * np.sum((2. * Llls + 1) * result[zi]['S'] * result[zi]['A0i_j'] * result[zi]['A0i_j'] / result[zi]['Aii'] / A00)
      
      for j, zj in enumerate(result):
        Fisher[i, j] += np.sum((2. * Llls + 1) * fsky * 2. * result[zi]['S'] ** 2. * result[zi]['r'] * result[zj]['r'] \
                                   * np.sqrt(1. / result[zi]['Aii'] / result[zj]['Aii']) * result[zi]['A0i_j'] * result[zj]['A0i_j'] / A00)  

    iFish    = np.linalg.inv(Fisher)
    diFish   = np.diag(iFish)

    for i, zz in enumerate(result):
        dstr  = "\tz:  %.2lf \t\t Schur-Limber:  %.2lf \t\t Limber:  %.2lf" % (zz, 100. * result[zz]['ferr_ii'],\
                                                                                   100. * np.sqrt(diFish[i])/result[zz]['Np'])

        dstr += "\t\t Schur at Lmin: %.2lf, and Lmax: %.2lf"                % (result[zi]['S'][0], result[zi]['S'][-1])

        if printit:
          print(dstr)

    ##  Define output. 
    output = []

    for outz in percentiles: 
        index  = np.where(np.abs((zs - outz)) == np.min(np.abs(zs - outz)))[0]

        ##  Save 
        output += [100. * result[zs[index][0]]['ferr_ii'], 100. * np.sqrt(diFish[index]) / result[zs[index][0]]['Np']]

    return  output


if __name__ == '__main__':
  import  os
  import  pylab              as      pl

  from    pickle             import  dump, load
  from    ilim               import  get_nbar_nocontam
  from    pz_tools           import  percentiles  
  from    schechter.gen_pz   import  peakz            as  _peakz
  from    schechter.get_shot import  get_shot
  from    schechter.get_pz   import  get_pz
  from    get_bz             import  bz_callmodel


  print('\n\nWelcome to a McQuinn and White clustering redshift forecaster.')

  evaluate     =   True
  add_desi     =  False

  dz           =   0.1
  
  fsky         =  0.01
  fover        =  1.00

  band         =   'r'

  setup        = {'BX':  {'colors': ['goldenrod', 'tan',         'y'], 'spec-maglim': 24.0, 'phot-maglim': 25.5,'decband': 'R', 'file': 'BXDrop.dat'},\
                   'u':  {'colors': ['darkblue',  'deepskyblue', 'b'], 'spec-maglim': 24.0, 'phot-maglim': 24.6,'decband': 'R', 'file': 'MalkanDrop.dat'},\
                   'g':  {'colors': ['darkgreen', 'limegreen',   'g'], 'spec-maglim': 25.5, 'phot-maglim': 25.8,'decband': 'i', 'file': 'gDrop.dat'},\
                   'r':  {'colors': ['darkred',   'indianred',   'r'], 'spec-maglim': 25.5, 'phot-maglim': 25.8,'decband': 'z', 'file': 'rDrop.dat'}}

  file         =  setup[band]['file']
  colors       =  setup[band]['colors']

  spec_mlim    =  setup[band]['spec-maglim']
  phot_mlim    =  setup[band]['phot-maglim']

  ps           =  get_pz(band)
  bs           =  lambda z:  bz_callmodel(z, spec_mlim)

  ##  Get the z percentiles for this dropout p(z).                                                                                                          
  percentiles = percentiles(ps, printit=True)

  pp           =  get_pz(band)

  peakz        =  _peakz(ps)

  zmin         =  peakz - 1.0
  zmax         =  peakz + 1.0

  ##  Get dropout Schechter photometric sample counts for given band. 
  root         =  os.environ['LBGCMB']
  data         =  np.loadtxt(root + '/dropouts/schechter/dat/%s' % file)

  ##  Number of photometric galaxies per sq. deg. for given detection band limit. 
  ms           =  data[:,0][::-1]
  Npz          =  data[:,1][::-1]

  ##  Cut to half mags. above 24.0                                                                                                                         
  valid        =  (ms % 0.5 == 0) & (ms >= 24.0) & (ms < 26.5)

  ms           =    ms[valid]
  Npz          =   Npz[valid]

  bps          =  [lambda z:  bz_callmodel(z, m) for m in ms]
  
  print('Evaluating for %s-dropout nphot:' % band)

  for i, m in enumerate(ms):
      print('%.2lf \t %.2lf \t %.2lf' % (m, bps[i](peakz), Npz[i]))

  print

  ##  S tends to infinite if Ns = Np in the shot noise limit.                                                                                                
  Nsz          =  np.logspace(2.0, np.log10(5000.), 5, base=10.)
  intlp_zs     =  []

  print('\nEvaluating for nspec:\n' + ''.join('%.2lf\n' % x for x in Nsz))

  if evaluate:
    ##  Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                                   
    cambx        =  CAMB()
    Pk_interps   =  get_PkInterps(cambx)

    ##  Input NLlls is ignored in the log10=False case.                                                                                                
    NLlls, Llls, nmodes = prep_Llls(NLlls = 60, Lmin = 60., Lmax = 5000., log10=False)

    results = []
     
    for Ns in Nsz:
      for ii, Np in enumerate(Npz):
          print("Solving for:  Ns: %.2lf \t mlim:  %.2lf \t Np: %.2lf" % (Ns, ms[ii], Np))

          result = Fisher(Pk_interps, Llls, Ns, Np, ps, bs, pp, bps[ii], dz=dz, fsky=fsky, zmin=zmin, zmax=zmax, fover=fover,\
                          percentiles=percentiles, intlp_zs=intlp_zs)

          results.append([Ns, ms[ii], Np] + result)

    results = np.array(results)

    np.savetxt(os.environ['LBGCMB'] + '/mqw/dat/mqw_%sdrops.txt' % band, results, fmt='%.4le', delimiter='\t')
    
  ##  And plot ...
  data = np.loadtxt(os.environ['LBGCMB'] + '/mqw/dat/mqw_%sdrops.txt' % band)
    
  Nsz  = np.unique(data[:,0])  ##  N spec. 
  mms  = np.unique(data[:,1])  ##  mag. lim.
  Npz  = np.unique(data[:,2])  ##  N phot.  

  ##  
  latexify(fig_height=2.2, columns=2, fontsize=12)
  
  if add_desi:
      fig, axs = plt.subplots(1, 4, sharey=True)
      index    = 1
  
  else:
      fig, axs = plt.subplots(1, 3, sharey=True)
      index    = 0 
      
  colors = ['k', 'b', 'r', 'indigo', 'y', 'sandybrown']

  for color, Np in zip(colors, Npz):
   dat = data[data[:,2] == Np]
  
   for kk, percentile in enumerate(percentiles):
     if  (kk == index) & (' (%.1lf)' % dat[0,1] in [' (24.0)', ' (24.5)']):
         label =  "%s" % (sci_notation(Np, decimal_digits=1, precision=None, exponent=np.int(np.floor(np.log10(Np))))) + ' (%.1lf)' % dat[0,1]    

     elif (kk == index + 1) & (' (%.1lf)' % dat[0,1] not in [' (24.0)', ' (24.5)']):
         label =  "%s" % (sci_notation(Np, decimal_digits=1, precision=None, exponent=np.int(np.floor(np.log10(Np))))) + ' (%.1lf)' % dat[0,1]

     else:
         label = ''

     axs[kk + index].semilogx(dat[:,0] / 1.e3, dat[:,3 + 2*kk], '-',  label = label, c=color, alpha=0.6, lw=1)

     ##  
     if   (kk == index) & (' (%.1lf)' % dat[0,1] in [' (24.0)', ' (24.5)']):
       axs[kk + index].legend(ncol=1, title=r'$z \simeq %.2lf$' % percentile, handlelength=.5, fontsize=7, title_fontsize=8)

     elif (kk == index + 1) & (' (%.1lf)' % dat[0,1] not in [' (24.0)', ' (24.5)']):
       axs[kk + index].legend(ncol=1, title=r'$z \simeq %.2lf$' % percentile, handlelength=.5, fontsize=7, title_fontsize=8) 

     else:
       axs[kk + index].legend(ncol=1, title=r'$z \simeq %.2lf$' % percentile, handlelength=.5, fontsize=7, title_fontsize=8)  

  for ax in axs:
    ax.fill_between(np.arange(0., 1.1e6, 1.e6), 0., 1., color=colors[0], alpha=0.2)
      
    ##  title  = r'$%.1lf < z < %.1lf$' % (zmin, zmax) + ' for ' + r'$f_{\rm{sky}} = %.2lf$, ' % fsky + 'd$z$=%.1lf' % dz 
    ##  title +=  ' and ' + r'$f_{\rm{over}} = %.1lf$' % fover
      
    ax.set_xlabel(r'$N_s \ [10^3]$', fontsize=8)
    ax.set_xscale('linear')

    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)

    ax.set_xlim(0.1, np.log10(5000.))
    ax.set_ylim(0.0, 8.000)

    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')

  axs[0].set_ylabel(r'$(\delta N_p \ / \ N_p) \ [\%]$', fontsize=8)

  if add_desi:
    axs[0].legend(ncol=1, title=r'$    (z \simeq %.2lf)$' % intlp_zs[0], handlelength=.5, fontsize=8)

  plt.tight_layout()

  pl.savefig(os.environ['LBGCMB'] + '/mqw/plots/mqw_%sdrops.pdf' % band)
  
  print('\n\nDone.\n\n')
