import  numpy  as  np


def prep_Llls(NLlls = 40, Lmin = 50., Lmax = 2000., log10=False):
  if log10:
    if Lmin < 25:
      ## Log spacing does not match nicely to integers for Lmin < 25.
      Llls    = np.concatenate([np.arange(Lmin, 25, 1), np.logspace(np.log10(25), np.log10(Lmax), NLlls - (25 - Lmin))])
    
    else:
      Llls    = np.logspace(np.log10(Lmin), np.log10(Lmax), NLlls)

  else:
    Llls      = np.arange(Lmin, Lmax + 1) 
    NLlls     = len(Llls)

  Llls        = np.ceil(Llls).astype(int)

  nmodes      = Llls - np.roll(Llls, 1)
  nmodes[0]   = 1

  return  NLlls, Llls, nmodes


if __name__ == '__main__':
  print('\n\nWelcome to prep. Llls.')

  NLlls, Llls, nmodes = prep_Llls(NLlls = 50, Lmin = 2., Lmax = 2000., log10=False)

  for i, L in enumerate(Llls):
    print L, nmodes[i], len(Llls), NLlls

  print('\n\nDone.\n\n')
