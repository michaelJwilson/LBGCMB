import  numpy        as      np 
import  pylab        as      pl

from    nbar         import  mlimitedM
from    reddy        import  samplestats
from    matplotlib   import  patches
from    schechterfn  import  SchechterMfn


stats       =  samplestats(printit = True)    ## Luminosity fn. of * all star forming galaxies *.                                                                                                                                         
                                              ## Note:  [\phi*] = [h_70/Mpc]^3 per mag, for M*_AB(1700 \AA).                                                                                                                             
alpha       =  -1.30                          ## Van Daalen and White. 
## alpha    =   0.50                          ## Convergent \bar for mlim -> oo 

M_star      =  -21.1
phi_star    =  stats['LBG']['phi_star']


def plot_figone():
  MM          =  np.linspace(M_star - 3., M_star + 6., 100)
  PhiMUV      =  SchechterMfn(MM, phi_star, M_star,  alpha)

  ##  Compare to Van Daalen and White, https://arxiv.org/pdf/1703.05326.pdf                                                                                                                                                              
  pl.plot(MM, np.log10(PhiMUV), 'k-')

  ##  Plot M_star.                                                                                                                                                                                                                       
  pl.axvline(M_star, ymin=0., ymax=1., c='k', lw=2, label=r'$M_{\star}$')

  for z, color in zip([0.1, 0.35, 0.60, 0.85, 1.10], ['b', 'c', 'g', 'y', 'r']):
    Mbr, Lbr = mlimitedM(z, 22.0, M_star)
    Mft, Lft = mlimitedM(z, 22.5, M_star)

    ax       = pl.gca()
    
    height   = np.log10(SchechterMfn((Mft + Mbr)/2., phi_star, M_star,  alpha))
    height   = 9.0 - height

    rect     = patches.Rectangle((Mft, -7.), width=(Mbr - Mft), height=height, linewidth=1, edgecolor=color, facecolor=color, alpha=0.4, label=str(z))

    ##  Add the patch to the Axes                                                                                                                                                                                                        
    ax.add_patch(rect)

  pl.xlim(-15., -23.)
  pl.ylim(-7.,   -2.)

  pl.xlabel(r'$M$')
  pl.ylabel(r'$\log_{10}|\Phi(M)|$')

  pl.title(r'Van Daalen & White, Fig. 1; No k-correction.')

  legend = pl.legend(ncol=2, loc=3)
  
  pl.savefig('plots/vandaalen_figone.pdf')


def plot_figtwo():
  from  cosmo   import  cosmo
  from  params  import  params
  from  nbar    import  dVols, comovdensity, projdensity


  dz     = 5e-2
  zs     = np.arange(0., 2.,    dz)
  
  zs, Vs = dVols(zs, cosmo, params)


  nzs    = []
  pnbars = []
  
  for zee in zs:
    nz    = comovdensity(zee + dz / 2., phi_star, M_star, alpha, type='app', mlim = 22.5, printit=False)
    nzs.append(10. ** nz)

    pnbar = projdensity(zee,  zee + dz, phi_star, M_star, alpha, mlim=22.5, printit = True)
    pnbars.append(pnbar)

  nzs    = np.array(nzs)
  pnbars = np.array(pnbars)

  dNdzs  = nzs * Vs 

  pl.plot(zs,  Vs/Vs.max(),        'k-',   label=r'$\Delta V(z) \ \ [(h^{-1} \rm{Mpc})^3$]')
  pl.plot(zs, nzs/nzs.max(),       'k--',  label=r'$\Delta \Phi(z)$')
  pl.plot(zs, dNdzs/dNdzs.max(),   'b-',   label=r'$dN/dz$')

  ## Projected density calc. for a thin dz slice = dN/dz * dz.  
  pl.plot(zs, pnbars/pnbars.max(), 'y--',  label=r'')

  pl.xlabel(r'$z$')

  pl.legend(loc=1)

  pl.savefig('plots/vandaalen_figtwo.pdf')

def plot_figthree():
  from  scipy.special  import  gamma
  from  nbar           import  dVols, comovdensity

  
  mags   = np.arange(18.0, 30., 0.1)  

  result = []

  for mag in mags:
    nbar = comovdensity(1.0, phi_star, M_star, alpha, type='app', mlim=mag, printit=False)

    result.append(nbar)
    
  result = np.array(result)
  pl.plot(mags, result)


  limit  = np.log10(phi_star * gamma(1. + alpha))
  pl.axhline(limit, c='k', linestyle='--', label= r'$\phi_{*} \ \Gamma(1 + \alpha)$')

  ## print(alpha, gamma(1. + alpha), phi_star * gamma(1. + alpha), limit) 

  pl.xlabel(r'$m_{\rm{lim}}$')
  pl.ylabel(r'$\log_{10} |\bar n|$')

  pl.ylim(-50., 0.)

  pl.legend(loc=4)

  pl.savefig('plots/vandaalen_figthree.pdf')


if __name__ == "__main__":
  print("\n\nReproducing Van Daalen & White sanity checks for dN/dz.\n\n")
  
  ## plot_figone()

  plot_figtwo()

  ## plot_figthree()

  print("\n\nDone.\n\n")
