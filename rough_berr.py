import  numpy  as  np


def rough_berr(Cgg, sigCgg, Ckg, sigCkg, samplevar_limit = False):
  """
  Given the variance of binned estimates of Ckg(z) and Ggg(z),                                                                                             
  roughly estimate the optimal error on b(z) (obtained by inverse variance weighting).                                                                     

  if     \hat b(ell) = \hat Cgg(ell) / \hat Ckg(ell) for each binned ell estimate.                                                                         
  then  (\sig b(ell) / b)^2  = (\sig Cgg(ell) / Cgg(ell))^2 + (\sig Ckg(ell) / Ckg(ell))^2  for each binned ell estimate.                                  

  Where inverse variance weighting of the binned ell estimates gives                                                                                  
         \var \hat b = 1. /  \sum_ell  1./ \var \hat b_l                                                                                                   
  """
  
  sig_bell_onb       = np.sqrt((sigCgg/Cgg)**2. + (sigCkg/Ckg)**2.)
  invvar_berr_onb2   = 1./np.sum(1./sig_bell_onb**2.)

  if samplevar_limit == False:
    print "With    shotnoise, fractional error on b is %.3lf percent" % (100.*np.sqrt(invvar_berr_onb2))

  else:
    print "Without shotnoise, fractional error on b is %.3lf percent" % (100.*np.sqrt(invvar_berr_onb2))

  return 100.*np.sqrt(invvar_berr_onb2)

def add_berrbox(frac_berrs):
  offsetbox  = TextArea(r"$z$       $\sigma_b/b$ [%]" + "".join("\n%.1lf    %.2lf" % (x, frac_berrs[x]) for x in frac_berrs.keys()) , minimumdescent=False)
  xy         = (0.675, 0.820)
  ab         = AnnotationBbox(offsetbox,
                              xy,
                              xybox=(0, 1),
                              xycoords='figure fraction',
                              boxcoords="offset points")
  ax         = plt.gca()

  ax.add_artist(ab)


if __name__ == "__main__":
    print("Welcome to rough estimate of the linear bias estimate.")


    print("Done.")
