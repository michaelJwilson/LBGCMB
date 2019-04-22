import  numpy   as      np
 
from    cosmo   import  cosmo
from    utils   import  comoving_distance
from    params  import  get_params


params = get_params()

if __name__ == '__main__':
    print('\n\nWelcome.\n\n')

    for zee in np.arange(0., 7., 1.):
        print('%.3lf \t %.3le \t %.3le' % (zee, comoving_distance(zee), params['h_100'] ** 3. * cosmo.comoving_volume(zee).value / 1.e9))

    print('\n\nDone.\n\n')
