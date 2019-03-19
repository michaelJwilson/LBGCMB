def densespec(Llls, zee, result, fsky=1.e-3, ellim=1.e3):
    ##  Eqn. (1) of McQuinn and White;  Fractional error on Np for bin i.                                                                                   
    result    = result[zee]['beta'][Llls < ellim]
    Llls      =                Llls[Llls < ellim]

    result    = 0.1 * (fsky / 1.e-3) ** -0.5 * (ellim / 1.e3) ** -1. * (result / 0.1)**-0.5
    norm      = np.sum(2. * Llls + 1.)

    return  np.sum((2. * Llls + 1.) * result) / norm

def sparsespec(Llls, zee, result, Nspec = 1.3, ellim=1.e3):
    ##  Eqn. (2) of McQuinn and White;  Fractional error on Np for bin i.                                                                                   
    result    = result[zee]['beta'][Llls < ellim]
    Llls      =                Llls[Llls < ellim]

    norm      = np.sum(2. * Llls + 1.)
    result    = np.sum((2. * Llls + 1.) * result) / norm

    return  (Nspec / 1.e3)** -0.5 * (result / 0.1)**-0.5
