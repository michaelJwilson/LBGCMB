from    pz2nbar            import  nbar_convert


def Ng(ilim, deg=True):
    ##  McQuinn and White, below eqn. (3).                                                                                                                  
    result  = 1.7 * 10. ** (5. + 0.31 * (ilim - 25.))    ##  [deg2]                                                                                          

    if deg:
      return  result

    else:
      return  nbar_convert(result, unit='str')           ##  [steradians2] 
