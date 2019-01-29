def get_params():
  ## Planck '18 + BAO.  
  ## params  = {'om_m': 0.3106, 'om_b': 0.04898, 'om_L': 0.6894, 'h_100': 0.6770, 'sig_8': 0.811, 'Tcmb': 2.7255, 'zscatter': 1100, 'ns': 0.96824, 'As': 2.1073e-9}

  ## Run PB (modi_lensing.py)
  params  = {'om_m': 0.292, 'om_b': 0.0468, 'om_L': 0.708, 'h_100': 0.69, 'sig_8': 0.819, 'Tcmb': 2.786, 'zscatter': 1100, 'ns': 0.9688, 'As': 2.187e-9}

  ## Dark sky cosmology:                                                                                                                                 
  ## params  = {'om_m':0.295, 'om_b':0.0468, 'om_L':0.705, 'h_100':0.688, 'sig_8':0.835, 'Tcmb':2.786, 'zscatter': 1100, 'ns': 0.9688, 'As': 2.187e-9}

  ## Add in cold dark matter as the difference between matter and baryons. neutrinos?
  params['om_cm'] = params['om_m'] - params['om_b'] 
  params['h_70']  = (10./7.)*params['h_100']

  return  params

def get_hbias():
  ## Parameters for Modi et al. cross-correlation fits as a function of redshift; Set limit of z-range of params to 4.                 
  hbias = { 1.0: {'b1':4.3801e-01, 'b2':2.4311e-01, 'bs':-2.6139e-01, 'ahm':-1.1136e+00, 'amm':-2.1587e-01}, \
            1.5: {'b1':9.8217e-01, 'b2':3.6117e-03, 'bs': 7.7140e-04, 'ahm': 6.8840e-01, 'amm': 2.9043e-03}, \
            2.0: {'b1':1.5872e+00, 'b2':1.1692e+00, 'bs': 1.4286e-01, 'ahm':-2.5484e+00, 'amm': 1.1870e+00}, \
            2.5: {'b1':2.3665e+00, 'b2':2.7414e+00, 'bs':-3.2667e+00, 'ahm':-1.6343e+00, 'amm':-9.7393e-01}, \
            3.0: {'b1':3.2109e+00, 'b2':6.6665e+00, 'bs': 1.1453e+00, 'ahm':-6.0501e+00, 'amm':-3.0775e+00}, \
            3.5: {'b1':0.0000e+00, 'b2':0.0000e+00, 'bs': 0.0000e+00, 'ahm': 0.00000+00, 'amm': 0.0000e+00}}

  return  hbias
