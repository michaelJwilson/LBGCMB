import  numpy             as  np
import  pylab             as  pl 
import  collections
import  astropy.constants as  const

from    rrtemplate_io     import read_template

def get_onine(mJy = False):
    dat = np.loadtxt("/Users/M.J.Wilson/work/LBG-CMB/dropouts/stellar_spectra/uko9v.dat")

    ls  = dat[:,0]
    Fl  = dat[:,1]

    vs  = (1e10/ls) * const.c.to('m/s').value   ## ls [AA]; vs [Hz]; c [m/s.]
    Fv  = ls * Fl / vs 

    return ls, Fl


if __name__ == "__main__":
    ls, Fl      = get_onine()

    Fl          = Fl[ls <= 1.0e4] 
    ls          = ls[ls <= 1.0e4]

    norm        = Fl[-1] 

    pl.semilogy(ls, Fl, label='Star-O')
    
    types       = collections.OrderedDict()

    ## Redrock stellar templates.
    keys        = ['Star-M',  'Star-K', 'Star-G', 'Star-F',  'Star-A', 'Star-B']

    for type in keys:
        i           = 0

        types[type] = read_template(type=type, printit=False)
        
        ls          = types[type]['wave']
        Fv          = types[type]['temp_%d' % i]                                  ## [1e-17 ergs/s/A/cm$^2]

        vs          = (1e10/ls) * const.c.to('m/s').value                         ## ls [AA]; vs [Hz]; c [m/s.]                                                
        Fl          = vs * Fv / ls

        Fl          = Fl[ls <= 1.0e4]
        ls          = ls[ls <= 1.0e4]

        Fl         *= norm / Fl[-1]

        pl.semilogy(ls, Fl, lw=0.2, label='%s' % type)
     
    pl.xlabel(r'Wavelength [$\AA$]')
    pl.ylabel(r'$F_{\lambda}$')
    
    pl.xlim(3.e3, 1.0e4)
    ## pl.ylim(0., 0.25)

    pl.legend(ncol=2)
    
    pl.savefig('./stellar_spectra.pdf')
