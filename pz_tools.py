import  numpy         as      np

from    whitebook_pz  import  whitebook_pz


def  standard_normal(z):
   ## Note:  not normalised. 
   return  np.exp(- z * z / 2.)

def  percentiles(pz, printit=False, asarray=True):
   ## Check normalisation. 
   dz   = 0.01

   zs   = np.arange(0.0, 10., dz)
   ps   = pz(zs)

   norm = np.sum(ps) * dz 
   ps  /= norm
   
   cs   = np.cumsum(ps) * dz

   Q1   = zs[cs > 0.25][0]
   Q2   = zs[cs > 0.50][0]
   Q3   = zs[cs > 0.75][0]

   if  printit:
      print  zs
      print  ps
      print  cs
   
   if asarray:
     return  np.array([Q1, Q2, Q3])

   else:
     return  [Q1, Q2, Q3] 


if __name__ == '__main__':
   print('\n\nWelcome to p(z) tools.\n\n') 

   zs = np.arange(0.0, 10., 1.0)

   pz = whitebook_pz(zs, ilim = 25.3)

   percentiles(standard_normal)
   ## percentiles(whitebook_pz)

   print('\n\nDone.\n\n')
