def Reddy_colourcut(zs, luv, u, g, r, i, z, y, type='BX', nocolourcut = False):
    if type == 'BM':
        ##  1.5 < z < 2.0                                                                                                                               
        crit  = (g - r) >= -0.2
        crit &= (u - g) >= (g - r)     - 0.1
        crit &= (g - r) <= 0.2*(u - g) + 0.4
        crit &= (u - g) <= (g - r)     + 0.2

    elif type == 'BX':
        ##  2.0 < z < 2.5                                                                                                                              
        crit  = (g - r) >= -0.2
        crit &= (u - g) >= (g - r)     + 0.2
        crit &= (g - r) <= 0.2*(u - g) + 0.4
        crit &= (u - g) <= (g - r)     + 1.0

    else:
        raise  ValueError("Specified Reddy colour selection is not available.")
