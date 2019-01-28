import numpy as np


def effective_depth(dropband='g', depth='W'):
    ##  Table 1 of https://arxiv.org/pdf/1704.06004.pdf.                                                                                                       
    ##  Note:  'Effective' depth in the dropout band.                                                                                                          
    if  dropband == 'g' and depth == 'W':
        return (26.43 + 26.35 + 26.38 + 26.39 + 26.47 + 26.31) / 6.0

    elif dropband == 'g' and depth == 'D':
        return (26.73 + 26.56 + 26.77 + 26.69) / 4.0

    elif dropband == 'g' and depth == 'UD':
        return (27.15 + 27.13) / 2.0

    elif dropband == 'r' and depth == 'W':
        return (25.93 + 25.88 + 25.95 + 25.96 + 26.04 + 25.87) / 6.0

    elif dropband == 'r' and depth == 'D':
        return (26.30 + 26.19 + 26.13 + 26.25) / 4.0

    else:
        raise ValueError('\n\nCombination of %s and %s is not available.' % (dropband, depth))
