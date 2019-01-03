"""
A dictionary containing simple parameters of current and future CMB experiments. 
"""
bolometers  = {'Planck': {'fsky': 1.0, 'thetab': 7.0, 'DeltaT': 30.0, 'iterative': False},\
               'AdvACT': {'fsky': 0.2, 'thetab': 1.5, 'DeltaT': 12.0, 'iterative': False},\
                   'SO': {'fsky': 0.4, 'thetab': 1.4, 'DeltaT':  7.0, 'iterative':  True},\
                'CMBS4': {'fsky': 0.4, 'thetab': 1.4, 'DeltaT':  1.0, 'iterative':  True},\
               ## 'CORE':   {'fsky': 0.8, 'thetab': 2.0, 'DeltaT':  2.0, 'iterative': False},\
               ## 
               ## Schmittfull & Seljak (2017);  https://arxiv.org/pdf/1710.09465.pdf
               ## 'SS17':   {'fsky': 0.5, 'thetab': 1.0, 'DeltaT':  1.0, 'iterative': True}
              }


if __name__ == '__main__':
    print('\n\nWelcome to bolometers.\n\n')

    dict = bolometers['SO']

    print('\n\nDone.\n\n')
