from    astropy.io import fits

import  numpy as np
import  glob
import  sys
import  os


def get_filenames(printit = True):
  for line in infiles:
    if not line.strip().startswith("#"):
      [name, root, path]   = line.split() 
         
      allfits[name]        = roots[root] + path
         
      if printit:
        print(allfits[name])

  return  allfits

def print_header(allfits, name, printhead = True, tile = '00098'):  
  fname                    = allfits[name]
  
  print("\nGetting:  %s\n" % fname)
  
  try:
    dat                    =  fits.open(fname)

  except:
    infiles    = open('infiles.txt', 'r')

    if name == "FIBERS":
      fname                = os.environ['CHALLENGE'] + '/fiberassign/tile_%s.fits' % tile
      dat                  = fits.open(fname)

    else:
      raise ValueError("\n\nCannot find: %s\n\n" % fname)

  hdrs                     = [x.header for x in dat]
  hdrkeys                  = list(dat[1].header.keys()) 

  if printhead:
    print(repr(dat[1].header))
    print('\n\n')
  
  for i, x in enumerate(hdrkeys):
    if "TTYPE" in x:
      hdr = dat[1].header[i]

      if "COEFF" not in hdr:
        print("".join(np.str(x).ljust(25) for x in [i, hdr] + [dat[1].data[i][hdr] for i in np.arange(len(dat[1].data))]))
  

if __name__ == "__main__":
  print("\n\nWelcome to seefits.\n\n")

  if '-i' in sys.argv:
    allfits         =  {}
    allfits['USER'] =  sys.argv[sys.argv.index('-i') + 1]
  
  else:
    roots           = {'DMOCKS': os.environ['DMOCKS'], 'CHALLENGE': os.environ['CHALLENGE'], 'REALWORLD': os.environ['REALWORLD']}
    allfits         = get_filenames(printit = False)

  print("Available keys:  ")

  for x in allfits.keys():
    print(x)
  
  ## for x in ["ELG", "TARGETSTRUTH", "TARGETS", "CONDITIONS", "FIBERS", "IN2OUT"]:
  for x in allfits.keys():
    print_header(allfits, x)

  print("\n\nDone.\n\n")
