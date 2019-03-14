import os
import pandas as pd


def get_pEBV(printit=False):
  names       = ['loEBV', 'hiEBV', 'p(z=2)', 'p(z=3)']

  root        = os.environ['LBGCMB']
  data        = pd.read_csv(root + '/dropouts/reddy/dat/tabfive.dat', sep='\s+', skiprows=0, names=names)
  data['EBV'] = 0.5 * (data['loEBV'].values + data['hiEBV'].values)

  if printit:
    for i, x in enumerate(data['EBV'].values):
      print('%.2lf \t %.4lf' % (x, data['p(z=3)'].values[i]))

  return  data

if __name__ == '__main__':
    get_pEBV(printit=True)
