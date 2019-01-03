
from __future__ import absolute_import, division


import argparse
import numpy as np
from numpy.linalg.linalg import LinAlgError
import astropy.io.fits as pyfits
from numpy.polynomial.legendre import legval,legfit

from desispec.io import read_image
from desiutil.log import get_logger
from desispec.trace_shifts import read_psf_and_traces,write_traces_in_psf,compute_dx_from_cross_dispersion_profiles,compute_dy_from_spectral_cross_correlation,monomials,polynomial_fit,compute_dy_using_boxcar_extraction,compute_dx_dy_using_psf,shift_ycoef_using_external_spectrum,recompute_legendre_coefficients



def parse(options=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Measures trace shifts from a preprocessed image and an input psf, and writes the modified trace coordinates in an output psf file to be used for extractions.
Two methods are implemented.
1) cross-correlation : dx shifts are measured from cross-dispersion profiles of traces. 
   dy shifts (wavelength calibration) are measured in two steps, an internal calibration determined from the cross-correlation of fiber spectra obtained from a resampled boxcar extraction with their median, and the final wavelength calibration is obtained from the  cross-correlation of the median fiber spectrum (after a second boxcar extraction) with an external spectrum (given with --spectrum option).
   This method is efficient for measuring trace shifts on science exposures with a large sky background.
2) forward model : dy,dy shifts are determined simultaneously by a forward modeling of the image around a given external list of lines (--lines option).
   This method is in principle statistically optimal, but it is slow and cannot be applied to blended and broad sky lines. It is useful to shift traces from arc lamp images (though specex does the same thing in C++).""")
    
    parser.add_argument('--image', type = str, default = None, required=True,
                        help = 'path of DESI preprocessed fits image')
    parser.add_argument('--psf', type = str, default = None, required=True,
                        help = 'path of DESI psf fits file')
    parser.add_argument('--lines', type = str, default = None, required=False,
                        help = 'path of lines ASCII file. Using this option changes the fit method.')
    parser.add_argument('--spectrum', type = str, default = None, required=False,
                        help = 'path to a spectrum ASCII file (e.g. the DESIMODEL sky spectrum)')
    parser.add_argument('--outpsf', type = str, default = None, required=True,
                        help = 'path of output PSF with shifted traces')
    parser.add_argument('--outoffsets', type = str, default = None, required=False,
                        help = 'path of output ASCII file with measured offsets for QA')
    parser.add_argument('--degxx', type = int, default = 2, required=False,
                        help = 'polynomial degree for x shifts along x')
    parser.add_argument('--degxy', type = int, default = 2, required=False,
                        help = 'polynomial degree for x shifts along y')
    parser.add_argument('--degyx', type = int, default = 2, required=False,
                        help = 'polynomial degree for y shifts along x')
    parser.add_argument('--degyy', type = int, default = 2, required=False,
                        help = 'polynomial degree for y shifts along y')
    parser.add_argument('--nfibers', type = int, default = None, required=False,
                        help = 'limit the number of fibers for debugging')
    
    args = None
    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args) :

    global psfs

    log=get_logger()

    log.info("starting")
    
    
    # read preprocessed image
    image=read_image(args.image)
    log.info("read image {}".format(args.image))
    if image.mask is not None :
        image.ivar *= (image.mask==0)

    xcoef=None
    ycoef=None
    psf=None
    wavemin=None
    wavemax=None
    nfibers=None
    lines=None

    psf,xcoef,ycoef,wavemin,wavemax = read_psf_and_traces(args.psf)
    nfibers=xcoef.shape[0]
    log.info("read PSF trace with xcoef.shape = {} , ycoef.shape = {} , and wavelength range {}:{}".format(xcoef.shape,ycoef.shape,int(wavemin),int(wavemax)))
    
    if args.lines is not None :
        log.info("We will fit the image using the psf model and lines")
        
        # read lines
        lines=np.loadtxt(args.lines,usecols=[0])
        ok=(lines>wavemin)&(lines<wavemax)
        log.info("read {} lines in {}, with {} of them in traces wavelength range".format(len(lines),args.lines,np.sum(ok)))
        lines=lines[ok]
        

    else :
        log.info("We will do an internal calibration of trace coordinates without using the psf shape in a first step")
        
        
        
    if args.nfibers is not None :
        nfibers = args.nfibers # FOR DEBUGGING

    fibers=np.arange(nfibers)

    if lines is not None :

        # use a forward modeling of the image
        # it's slower and works only for individual lines
        # it's in principle more accurate
        # but gives systematic residuals for complex spectra like the sky
        
        x,y,dx,ex,dy,ey,fiber_xy,wave_xy=compute_dx_dy_using_psf(psf,image,fibers,lines)
        x_for_dx=x
        y_for_dx=y
        fiber_for_dx=fiber_xy
        wave_for_dx=wave_xy
        x_for_dy=x
        y_for_dy=y
        fiber_for_dy=fiber_xy
        wave_for_dy=wave_xy
        
    else :
        
        # internal calibration method that does not use the psf
        # nor a prior set of lines. this method is much faster
        
        # measure x shifts
        x_for_dx,y_for_dx,dx,ex,fiber_for_dx,wave_for_dx = compute_dx_from_cross_dispersion_profiles(xcoef,ycoef,wavemin,wavemax, image=image, fibers=fibers, width=7, deg=args.degxy)
        # measure y shifts
        x_for_dy,y_for_dy,dy,ey,fiber_for_dy,wave_for_dy = compute_dy_using_boxcar_extraction(xcoef,ycoef,wavemin,wavemax, image=image, fibers=fibers, width=7)

    
    # write this for debugging
    if args.outoffsets :
        file=open(args.outoffsets,"w")
        file.write("# axis f x y delta error\n")
        for e in range(dy.size) :
            file.write("0 %f %d %f %f %f %f\n"%(wave_for_dy[e],fiber_for_dy[e],x_for_dy[e],y_for_dy[e],dy[e],ey[e]))
        for e in range(dx.size) :
            file.write("1 %f %d %f %f %f %f\n"%(wave_for_dx[e],fiber_for_dx[e],x_for_dx[e],y_for_dx[e],dx[e],ex[e]))
            
        file.close()
        log.info("wrote offsets in ASCII file %s"%args.outoffsets)
    
    
    
    
    degxx=args.degxx
    degxy=args.degxy
    degyx=args.degyx
    degyy=args.degyy
    
    while(True) : # loop because polynomial degrees could be reduced
        
        log.info("polynomial fit of measured offsets with degx=(%d,%d) degy=(%d,%d)"%(degxx,degxy,degyx,degyy))
        try :
            dx_coeff,dx_coeff_covariance,dx_errorfloor,dx_mod,dx_mask=polynomial_fit(z=dx,ez=ex,xx=x_for_dx,yy=y_for_dx,degx=degxx,degy=degxy)
            dy_coeff,dy_coeff_covariance,dy_errorfloor,dy_mod,dy_mask=polynomial_fit(z=dy,ez=ey,xx=x_for_dy,yy=y_for_dy,degx=degyx,degy=degyy)
            
            log.info("dx dy error floor = %4.3f %4.3f pixels"%(dx_errorfloor,dy_errorfloor))

            log.info("check fit uncertainties are ok on edge of CCD")
            
            merr=0.
            for fiber in [0,nfibers-1] :
                for rw in [-1,1] :
                    tx = legval(rw,xcoef[fiber])
                    ty = legval(rw,ycoef[fiber])
                    m=monomials(tx,ty,degxx,degxy)
                    tdx=np.inner(dx_coeff,m)
                    tsx=np.sqrt(np.inner(m,dx_coeff_covariance.dot(m)))
                    m=monomials(tx,ty,degyx,degyy)
                    tdy=np.inner(dy_coeff,m)
                    tsy=np.sqrt(np.inner(m,dy_coeff_covariance.dot(m)))
                    merr=max(merr,tdx)
                    merr=max(merr,tdy)
                    #log.info("fiber=%d wave=%dA x=%d y=%d dx=%4.3f+-%4.3f dy=%4.3f+-%4.3f"%(fiber,int(wave),int(x),int(y),dx,sx,dy,sy))
            log.info("max edge shift error = %4.3f pixels"%merr)
            if degxx==0 and degxy==0 and degyx==0 and degyy==0 :
                break
        
        except ( LinAlgError , ValueError ) :
            log.warning("polynomial fit failed with degx=(%d,%d) degy=(%d,%d)"%(degxx,degxy,degyx,degyy))
            if degxx==0 and degxy==0 and degyx==0 and degyy==0 :
                log.error("polynomial degrees are already 0. we can fit the offsets")
                raise RuntimeError("polynomial degrees are already 0. we can fit the offsets")
            merr = 100000. # this will lower the pol. degree.
        
        if merr > 0.05 :
            if merr != 100000. :
                log.warning("max edge shift error = %4.3f pixels is too large, reducing degrees"%merr)
            
            if degxy>0 and degyy>0 and degxy>degxx and degyy>degyx : # first along wavelength
                if degxy>0 : degxy-=1
                if degxy>0 : degyy-=1
            else : # then along fiber
                if degxx>0 : degxx-=1
                if degyx>0 : degyx-=1
        else :
            # error is ok, so we quit the loop
            break
    
    # print central shift
    mx=np.median(x_for_dx)
    my=np.median(y_for_dx)
    m=monomials(mx,my,degxx,degxy)
    mdx=np.inner(dx_coeff,m)
    mex=np.sqrt(np.inner(m,dx_coeff_covariance.dot(m)))

    mx=np.median(x_for_dy)
    my=np.median(y_for_dy)
    m=monomials(mx,my,degyx,degyy)
    mdy=np.inner(dy_coeff,m)
    mey=np.sqrt(np.inner(m,dy_coeff_covariance.dot(m)))    
    log.info("central shifts dx = %4.3f +- %4.3f dy = %4.3f +- %4.3f "%(mdx,mex,mdy,mey))
    
    # for each fiber, apply offsets and recompute legendre polynomial
    log.info("for each fiber, apply offsets and recompute legendre polynomial")
    xcoef,ycoef = recompute_legendre_coefficients(xcoef=xcoef,ycoef=ycoef,wavemin=wavemin,wavemax=wavemax,degxx=degxx,degxy=degxy,degyx=degyx,degyy=degyy,dx_coeff=dx_coeff,dy_coeff=dy_coeff)
    
    
    # use an input spectrum as an external calibration of wavelength
    if args.spectrum :
         
        
        log.info("write and reread PSF to be sure predetermined shifts were propagated")
        write_traces_in_psf(args.psf,args.outpsf,xcoef,ycoef,wavemin,wavemax)
        psf,xcoef,ycoef,wavemin,wavemax = read_psf_and_traces(args.outpsf)
                
        ycoef=shift_ycoef_using_external_spectrum(psf=psf,xcoef=xcoef,ycoef=ycoef,wavemin=wavemin,wavemax=wavemax,
                                                  image=image,fibers=fibers,spectrum_filename=args.spectrum,degyy=args.degyy,width=7)
    
    
        write_traces_in_psf(args.psf,args.outpsf,xcoef,ycoef,wavemin,wavemax)
        log.info("wrote modified PSF in %s"%args.outpsf)
        
    else :
        
        write_traces_in_psf(args.psf,args.outpsf,xcoef,ycoef,wavemin,wavemax)
        log.info("wrote modified PSF in %s"%args.outpsf)
        
    
