import numpy as np
import pylab as pl


print("\n\nWelcome to air2vacuum_lamdas.\n\n")

"""
Taken from http://classic.sdss.org/dr7/products/spectra/vacwavelength.html

Line     Air       Vacuum
H-beta   4861.363  4862.721
[O III]  4958.911  4960.295
[O III]  5006.843  5008.239
[N II]   6548.05   6549.86
H-alpha  6562.801  6564.614
[N II]   6583.45   6585.27
[S II]   6716.44   6718.29
[S II]   6730.82   6732.68
"""

VAC = 4862.721  

AIR = VAC / (1.0 + 2.735182e-4 + 131.4182 / VAC**2. + 2.76249e8 / VAC ** 4.)


print("\n\nVacuum:  %.6lf" % VAC)
print("\n\nAir:     %.6lf" % AIR)

print("\n\nDone.\n\n")
