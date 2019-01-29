import  os
import  numpy           as      np
import  astropy.units   as      u
import  pylab           as      pl

from    utils           import  latexify
from    prep_filters    import  prep_filters
from    BC03_maker      import  BC03_maker
from    app_mags        import  get_appmags
from    sklearn         import  linear_model
from    sklearn.metrics import  mean_squared_error, r2_score

##  latexify(fig_width=None, fig_height=None, columns=2, equal=False)                                                                                                                                                                     

print('\n\nWelcome to (BC03) LBG maker.')

ngal              =     2000  ##  Over written if rest-frame.                                                                                                                                                                              
dband             =      'r'
test              =   False
target_type       =   'BC03'
save              =   False
restframe         =   False

redshifts         =   3.0 + np.linspace(0.,   2.0,    20)
magnitudes        =  20.0 +   np.arange(0.,   7.0,   0.1)

flux, wave, meta  =  BC03_maker(ngal=ngal, restframe=restframe, printit=False, test=test, redshifts=redshifts, magnitudes=magnitudes, dband=dband, alliseds=True, calzetti=True, madau=True)

uwave             =  wave * u.AA
vs                =  uwave.to(u.Hz, equivalencies=u.spectral())

print
print(wave.shape)
print(flux.shape)
print(meta)


filters     =  prep_filters(['LSST', 'STEIDEL', 'SUBARU', 'JKC', 'HUBBLE'])

wave        =  wave * u.AA
vs          =  wave.to(u.Hz, equivalencies=u.spectral())

## y=x.
xs          = np.arange(-1.e2, 1.e2, 0.1)

pl.plot(xs, xs, 'k-')  

X = []

for i, x in enumerate(flux):  
  interim   =  x * u.erg / u.s / u.cm / u.cm / u.AA
  Fv        =  interim.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies = u.spectral_density(wave))
  Fv       *=  1.e-17

  mags      =  get_appmags(vs.value, Fv.value, filters, printit = False)

  ACS_beta  =  5.30 * (mags['ACS_F775W'] - mags['ACS_F850LP']) - 2.04
  LSST_beta =  5.30 * (mags['i'] - mags['z']) - 2.04

  X.append([ACS_beta, LSST_beta])

  pl.plot(ACS_beta, LSST_beta, 'ko', markersize=3)

X           = np.array(X)

X1          = X[:,0].reshape(-1, 1)
X2          = X1 ** 2.

Y1          = X[:,1].reshape(-1, 1)

##  Linear regression fit.
reg         = linear_model.LinearRegression(fit_intercept=True).fit(X1, Y1)

print(reg.score(X1, Y1))
print(reg.coef_)
print(reg.intercept_)

rms_res     = np.std(Y1 - reg.predict(X1))

print(rms_res)

xs          = np.arange(-5., 10., 0.1).reshape(-1, 1)
pl.plot(xs, reg.predict(xs), 'r-', label='Y = %.3lf * X + %.3lf (rms. residual:  %.3lf)' % (reg.coef_[0][0], reg.intercept_[0], rms_res))

pl.xlim(-5., 10.)
pl.ylim(-5., 10.)

pl.xlabel(r'ACS  $\beta$')
pl.ylabel(r'LSST $\beta$')
pl.legend()

##  pl.show()
pl.savefig(os.environ['BEAST'] + '/gal_maker/plots/ACS-beta.pdf')
