from   astropy.cosmology  import FlatLambdaCDM
from   params             import get_params


params = get_params()
cosmo  = FlatLambdaCDM(H0=100.*params['h_100'], Om0=params['om_m'], Ob0=params['om_b'])

if __name__ == '__main__':
    print('\n\nWelcome to a cosmo instance generator.\n\n')

    print('\n\nDone.\n\n')
