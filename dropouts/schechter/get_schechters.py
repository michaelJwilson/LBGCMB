def get_schechters(stats, key):
    '''
    Given a dropout stats dictionary,
    and a key to the dropout type,
    return z, alpha, M*, Phi*. 
    '''

    if key == 'u':
        key = 'Malkan'
    
    midz     =  stats[key]['z']
    alpha    =  stats[key]['schechter']['alpha']
    M_star   =  stats[key]['schechter']['M_star']
    phi_star =  stats[key]['schechter']['phi_star']

    return  midz, alpha, M_star, phi_star 


if __name__ == '__main__':
    from    reddy         import  samplestats as rsamplestats
    from    specs         import  samplestats as gsample_stats
    from    Malkan.specs  import  samplestats as usample_stats
    

    print('\n\nWelcome to get_schechters.\n\n')

    key   = 'LBG'
    stats =  rsamplestats(printit = True)

    midz, alpha, M_star, phi_star = get_schechters(stats, key)

    print(midz, alpha, M_star, phi_star)

    print('\n\nDone.\n\n')
