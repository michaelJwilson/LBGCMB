import  os
import  numpy              as  np
import  pylab              as  pl
import  matplotlib.pyplot  as  plt
import  threedhst.eazyPy   as  eazy
import  threedhst.utils    as  utils

from    astropy.table      import  Table


def write_table(results, header, depths, blockcount, ddir='scratch'):
  t = Table(rows=results, names=tuple(header))

  if depths is 'Full':
    t.write('colors/%s/3DHST_%s_Full_%d.fits' % (ddir, field, blockcount), format='fits', overwrite=True)

  elif depths is None:
    t.write('colors/%s/3DHST_%s_DefaultDepths_%d.fits' % (ddir, field, blockcount), format='fits', overwrite=True)

  else:
    t.write('colors/%s/3DHST_%s_CustomDepth_%d.fits' % (ddir, field, blockcount), format='fits', overwrite=True)


if __name__ == '__main__':
  ddir        = 'scratch'

  fnu         =   True
  save        =  False
  plotit      =   True

  ##  fields  = ['aegis',  'cosmos', 'uds', 'goodsn'] ## ['goodss']  ## master catalogue also available.
  fields      = ['UVUDF']

  ngal        =  9000
  blocksize   =   500       ##  Partition output, e.g. 9000, 500.

  ##  'Five sigma depth' for degraded imaging; noise realisation of uniform noise in Fv.
  for depths in ['Full']:   ## ['Full', None]. 
    for field in fields:
      results    = []
      blockcount =  1

      ## `idx` is zero-indexed.
      for idx in np.arange(ngal):
        if field == 'UVUDF':
          EAZYROOT         = '/global/homes/m/mjwilson/Rafelski15/UVUDF-EAZY/'
          MAIN_OUTPUT_FILE = 'photz'

        elif field == 'uds':
          EAZYROOT         = '/global/homes/m/mjwilson/3D-HST/%s.v4.2/Eazy/' % field
          MAIN_OUTPUT_FILE = '%s_3dhst.v4.2' % field
        
        else:          
          EAZYROOT         = '/global/homes/m/mjwilson/3D-HST/%s.v4.1/Eazy/' % field
          MAIN_OUTPUT_FILE = '%s_3dhst.v4.1' % field


        header, result = eazy.plotExampleSED(idx=idx, pdfname='plot.pdf', MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE,\
                                             OUTPUT_DIRECTORY=EAZYROOT, CACHE_FILE='Same', lrange=[1.e3, 1.e5], axes=None, individual_templates=False,\
                                             fnu=fnu, depths=depths, field=field, plotit=plotit)

        results.append(tuple(result))
        
        if save & (idx % blocksize == 0):
          write_table(results, header, depths, blockcount, ddir=ddir)

          ##  Wipe results.
          results     = []
          blockcount +=  1

      if save:
        ##  Catch the rest. 
        write_table(results, header, depths, blockcount, ddir=ddir)

      else:
        print('  '.join(item.ljust(6) for item in header))
        print(np.array(results))

print("\n\nDone.\n\n")
