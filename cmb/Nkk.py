import  os
import  numpy              as      np
import  matplotlib.pyplot  as      plt

from    prep_camb          import  Clxy
from    scipy.integrate    import  simps
from    pickle             import  load, dump
from    pmh                import  Pmm, get_PkInterps


plt.style.use('ggplot')

def prep_lmaxmask(nu, nv):
  mask_3000 = np.ones_like(nu)
  mask_5000 = np.ones_like(nv)

  mask_3000[nu > 3000] = 0.0
  mask_3000[nv > 3000] = 0.0

  mask_5000[nu > 5000] = 0.0
  mask_5000[nv > 5000] = 0.0

  mask_3000[nv <    2] = 0.0
  mask_5000[nv <    2] = 0.0

  return {'TT': mask_3000, 'TE': mask_3000, 'EE': mask_5000, 'BB': mask_5000, 'TB': mask_5000, 'EB': mask_5000}

def outer(*vs):  ## Form the outer product of N 1-dim arrays. 
    return reduce(np.multiply.outer, vs)

def Nkk(lensCl_interps, nolensCl_interps, Ls, terms=['TT', 'TE', 'EE', 'EB'], thetab=1., DeltaT=1., iterative=False, pickleit=True):
  '''
  Ckk noise curve for the flat-sky limit, appropriate to high ell. 
  see  Hu and Okamoto 2002, arXiv 0111606
  '''
  try:
    root    = os.environ['LBGCMB']

    fpath   = root + '/pickle/Nkk_thetab_%.2lf_DeltaT_%.2lf_' % (thetab, DeltaT) + 'Nlls_%d_Lmin_%d_Lmax_%d_' % (len(Ls), Ls[0], Ls[-1])
    fpath  += '_'.join(terms) + '_iterative_%d.p' % np.int(iterative == True)

    result  = load(open(fpath, 'rb'))

    print('\n\nLoaded: %s successfully.\n\n' % fpath)

    return  result

  except:
    print("Failed to load Nkk from %s, calculating from scratch.  Creating:\n" % fpath)
    print('%d < L < %d, in %d bins.' % (Ls[0], Ls[-1], len(Ls)))
    print('For:  thetab = %.3lf, DeltaT = %.3lf and with iteration = %d' % (thetab, DeltaT, np.int(iterative == True)))
    print('and:  ' + '_'.join(terms))

    lmax       =  5000.  ## NOTE: Needs changed depending on alpha.                                                                                           
    Nphi       =   200
    Nell       = 10000

    dphi       = 2.*np.pi/Nphi
    dell       =     lmax/Nell

    ell        = np.arange(1.0,     lmax + dell, dell)
    pphi       = np.arange(0.0, 2.*np.pi + dphi, dphi)

    l, phi, L  = np.meshgrid(ell, pphi, Ls)
  
    ux         = l*np.cos(phi)
    uy         = l*np.sin(phi)

    vx         = L - ux                       # Assumes \vec L = L * \hat x; L = u + v                                                                      
    vy         =    -uy
  
    phi_u      = np.arctan2(uy, ux)           # arctan2(y, x) returns [-pi, pi)                                                                               
    phi_v      = np.arctan2(vy, vx)           # Assumes vectors.

    phi_u[phi_u < 0.] += 2.*np.pi             # returns [0, 2 * pi]                                                                                          
    phi_v[phi_v < 0.] += 2.*np.pi

    phi_uv     = phi_u - phi_v                # \phi_{u, v} = \phi_u - \phi_v   

    nu         = (ux**2. + uy**2.)**0.5       # norm of u   
    nv         = (vx**2. + vy**2.)**0.5       # norm of v
  
    Lx         = ux + vx                      # \vec L   = \vec u + \vec v; by construction will lie along the x axis. 
  
    dLu        = ux*Lx                        # dot(L, u)
    dLv        = vx*Lx                        # dot(L, v)

    mask       = prep_lmaxmask(nu, nv)

    result     = np.zeros_like(ux)
  
    # All fa terms are unlensed;  TT, EE, BB have an added two in the denominator.
    if 'TT' in terms:
      print("Calculating TT.")

      fa = nolensCl_interps['TT'](nu)*dLu + nolensCl_interps['TT'](nv)*dLv    
      Fa = fa/(2. * Clxy(lensCl_interps, nu, 'TT', thetab, DeltaT) * Clxy(lensCl_interps, nv, 'TT', thetab, DeltaT))
     
      result += fa*Fa*mask['TT']
  
    if 'TE' in terms:                # cos phi_u = \hat x \cdot \hat u  
      print("Calculating TE.")

      c2phi_uv = np.cos(2.*phi_uv)

      # lensed                                                                                                                                           
      Cl_TT1  = Clxy(lensCl_interps, nu, 'TT', thetab, DeltaT)
      Cl_TT2  = Clxy(lensCl_interps, nv, 'TT', thetab, DeltaT)

      Cl_EE1  = Clxy(lensCl_interps, nu, 'EE', thetab, DeltaT)
      Cl_EE2  = Clxy(lensCl_interps, nv, 'EE', thetab, DeltaT)

      Cl_TE1  = Clxy(lensCl_interps, nu, 'TE', thetab, DeltaT)
      Cl_TE2  = Clxy(lensCl_interps, nv, 'TE', thetab, DeltaT)
                                                                                                # Note v, u ordering!                                       
      fa_uv   = nolensCl_interps['TE'](nu)*c2phi_uv*dLu + nolensCl_interps['TE'](nv)*dLv
      fa_vu   = nolensCl_interps['TE'](nu)*dLu          + nolensCl_interps['TE'](nv)*c2phi_uv*dLv

                                                                                               # Note v, u ordering! 
      num     = Cl_EE1*Cl_TT2*fa_uv - Cl_TE1*Cl_TE2*fa_vu
      denom   = Cl_TT1*Cl_EE2*Cl_EE1*Cl_TT2 - (Cl_TE1*Cl_TE2)**2.

      Fa      = num/denom
         
      result += fa_uv*Fa*mask['TE']
    
    if 'EE' in terms:
      print("Calculating EE.")

      fa = (nolensCl_interps['EE'](nu)*dLu + nolensCl_interps['EE'](nv)*dLv)*np.cos(2.*phi_uv)
      Fa = fa/(2. * Clxy(lensCl_interps, nu, 'EE', thetab, DeltaT) * Clxy(lensCl_interps, nv, 'EE', thetab, DeltaT))
      
      result += fa*Fa*mask['EE']

    if 'EB' in terms:
      print("Calculating EB.")

      fa = (nolensCl_interps['EE'](nu)*dLu - nolensCl_interps['BB'](nv)*dLv)*np.sin(2.*phi_uv)
      Fa = fa/(Clxy(lensCl_interps, nu, 'EE', thetab, DeltaT) * Clxy(lensCl_interps, nv, 'BB', thetab, DeltaT))

      if iterative:
        ## Iterative reconstruction accounts for a factor 2.5 for EB; Schmittfull & Seljak (2017).
        result += fa*Fa*mask['EB'] * 2.5

      else:
        result += fa*Fa*mask['EB']

    if 'BB' in terms:
      print("Calculating BB.")

      fa =  (nolensCl_interps['BB'](nu)*dLu + nolensCl_interps['BB'](nv)*dLv)*np.cos(2.*phi_uv)
      Fa = fa/(2. * Clxy(lensCl_interps, nu, 'BB', thetab, DeltaT) * Clxy(lensCl_interps, nv, 'BB', thetab, DeltaT))

      result += fa*Fa*mask['BB']
 
    if 'TB' in terms:
      print("Calculating TB.")

      fa = nolensCl_interps['TB'](nu)*np.sin(2.*phi_uv)*dLu
      Fa = fa/(Clxy(lensCl_interps, nu, 'TT', thetab, DeltaT) * Clxy(lensCl_interps, nv, 'BB', thetab, DeltaT))

      result  += fa*Fa*mask['TB']
 
    result    *= l

    result     = simps(result,  dx=dell, axis=0)
    result     = simps(result,  dx=dphi, axis=0)

    result    /= (2.*np.pi)**2.
    result     = 1./result

    result    *= (Ls*(Ls + 1.)/2.)**2.
  
    if pickleit:
      root     = os.environ['LBGCMB']

      fpath    = root + '/pickle/Nkk_thetab_%.2lf_DeltaT_%.2lf_' % (thetab, DeltaT) + 'Nlls_%d_Lmin_%d_Lmax_%d_' % (len(Ls), Ls[0], Ls[-1])
      fpath   += '_'.join(terms) + '_iterative_%d.p' % np.int(iterative == True)

      print('\n\nWriting: %s.' % fpath)

      dump(result, open(fpath, 'wb'))

    return  result
  

if __name__ == '__main__':
  import  pylab       as      pl 

  from    prep_Llls   import  prep_Llls
  from    prep_camb   import  CAMB
  from    utils       import  prefactor
  from    bolometers  import  bolometers  
  from    lensing     import  Ckk, var_Ckk
  from    snr         import  snr


  print('\n\nWelcome to Nkk.\n\n')

  ## Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                                         
  cambx                =  CAMB()

  Pk_interps           =  get_PkInterps(cambx)

  NLlls, Llls, nmodes  =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)

  ## No Detector noise -- this should be handled by Clxy of prep_camb.                                                                                     
  (lensCl_interps, nolensCl_interps) = cambx.get_Cls()

  ckk                  = Ckk(Pk_interps, Llls)

  pl.loglog(Llls, ckk, 'r')

  ## cmbexp            = 'SS17'                                                                                                                            
  for cmbexp in ['CMBS4']:
    fsky, thetab, DeltaT, iterative = bolometers[cmbexp]['fsky'], bolometers[cmbexp]['thetab'],\
                                      bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

    print('Calculating Nkk for %s' % cmbexp)
    
    nkk                = Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'], thetab=thetab,\
                             DeltaT=DeltaT, iterative=iterative, pickleit=True)

    ##  vkk            = var_Ckk(Llls, fsky, nkk, Pk_interps, samplevar_lim=False)

    ##  print('\n\nTotal S/N:  %.3lf' % snr(ckk, vkk, nmodes))

    pl.loglog(Llls, nkk, 'k', label=cmbexp)
    
  pl.xlim(2.,    2.e+3)
  pl.ylim(1.e-9, 4.e-6)

  pl.xlabel(r'$L$')
  ## pl.ylabel(r'$[ L (L + 1) ]^2 \ N_{\kappa \kappa} \ / \ (2 \pi) \ \ [ \times 10^{7} ]$')

  pl.savefig('plots/Nkk.pdf', bbox_inches='tight')

  print('\n\nDone.\n\n')
