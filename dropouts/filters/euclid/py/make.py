import numpy as np
import pylab as pl


##  Fig. 2 & Table 2 of https://arxiv.org/pdf/1710.09585.pdf

wave    = np.arange(4000., 22000., 1.0)
filters = {'VIS': [7150., 3550.], 'Y': [10850., 2750.], 'J': [13750., 4300.], 'H': [17725., 5250.]}

for filter in filters:
    [mid, fwid] = filters[filter]

    hwid        = fwid / 2.
    Ts          = np.zeros_like(wave)

    Ts[wave > (mid - hwid)] = 1.0

    if filter == 'VIS':
      Ts[wave > 7.e3]      *= 1. + 0.02 * (1. - (wave[wave > 7.e3] / 7.e3) ** 10.0)

    if filter == 'Y':
        x                   = 3. * 0.5 * np.pi * (wave - 8900.) / 2500.
        Ts                 += 0.025 * (np.sin(x) - 1.0)

    if filter == 'J':
        Ts                 *= np.exp(-((wave - mid) / fwid / 3.)**2.)

    if filter == 'H':
        Ts                 *= np.exp(-((wave - mid - 1000.) / fwid / 3.)**2.)

    Ts[wave < (mid - hwid)] = 0.0
    Ts[wave > (mid + hwid)] = 0.0

    cut = (mid - 1.05 * hwid < wave) & (wave < mid + 1.05 * hwid)

    pl.plot(wave[cut], Ts[cut], label=filter)

    np.savetxt('../%s.pb' % filter, np.c_[wave[cut], Ts[cut]], header='Fig. 2 & Table 2 of https://arxiv.org/pdf/1710.09585.pdf', fmt='%.6le')

pl.ylim(0.0, 1.2)
pl.xlabel(r'Wavelength [$\AA$]')
pl.ylabel(r'$T(\lambda)$')
pl.legend(loc=2, ncol=4, frameon=False)
pl.show()
