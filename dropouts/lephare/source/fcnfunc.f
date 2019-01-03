cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     function call by minuit to compute the chi2 between predicted 
c     and observed apparent magnitudes
      subroutine fcnfunc(npar,gin,f,x,iflag)

      implicit none

      integer*4 nbf,zadapt,bdincl,imagm
      INCLUDE 'dim_filt.decl'
      parameter (zadapt=30000)


      integer*4 i,k,ngals_ada,iflag,npar,fl1,fl2,fl_auto,nused,meth_ada
      real*8 f,a0,a1,a2,a3,sig,gin,x(npar)
      real*8 ab_ada(zadapt,nbf),sab_ada(zadapt,nbf)
      real*8 magm_ada(zadapt,nbf),zs_ada(zadapt)
c      integer*4 cont_ada(zadapt),mod_ada(zadapt)
      integer*4 mod_ada(zadapt)
      real*8 cont_ada(zadapt)
      real*8 chi,new_mag,po
      real*8 flux1,flux2,err_flux1

      external bdincl 

      common/func_int/k,ngals_ada,fl1,fl2,fl_auto,meth_ada,imagm
      common/func_real/sig
      common/func_tab/ab_ada,sab_ada,magm_ada,zs_ada
      common/func_tab2/cont_ada,mod_ada


      nused=0

      a0=x(1)
      a1=x(2)
      a2=x(3)
      a3=x(4)

c     initialise chi2
      chi=0.d0

c     loop over galaxies choose to adapt library
      do i=1,ngals_ada 

c      Check if we can use the apparent magnitude in this band
       if(bdincl(k-1,cont_ada(i),imagm-1).eq.1.or.
     .    cont_ada(i).eq.0)then

        po=0.d0
cccccccccccccccccccccccccccccccccccccccccccccccc
        if (meth_ada.eq.1) then 
c       Function of predicted color
          if(dabs(magm_ada(i,fl1)).le.50.d0.and.
     .      dabs(magm_ada(i,fl2)).le.50.d0) then
              po=magm_ada(i,fl1)-magm_ada(i,fl2)
          endif
        elseif (meth_ada.eq.2) then 
c       Redshift
          po=zs_ada(i)
        elseif (meth_ada.eq.3) then        
c       Modele
          po=dble(mod_ada(i))
        endif

c       New magnitudes with this coeficients 
        if(dabs(magm_ada(i,k)).le.50.d0)then
           new_mag=magm_ada(i,k)+a0+a1*po+a2*po**2.+a3*po**3.  
        else 
           new_mag=999.9
        endif


c       3 sigma clipping
        if(dabs(magm_ada(i,k)-ab_ada(i,k)).ge.3.*sig)new_mag=999.9

  
c       compute chi2
        if(dabs(ab_ada(i,k)).le.50.d0.and.sab_ada(i,k).le.0.6d0.and.
     .      sab_ada(i,k).gt.0.d0.and.new_mag.le.50.d0)then
c in mag 
c            chi=chi+((ab_ada(i,k)-new_mag)/sab_ada(i,k))**2.
c in flux
            flux1=10**(-0.4*(ab_ada(i,k)+48.59))
            flux2=10**(-0.4*(new_mag+48.59))
            err_flux1=flux1*sab_ada(i,k)/1.086
            chi=chi+((flux1-flux2)/err_flux1)**2.

            nused=nused+1
        endif

       endif

      enddo

      if(nused.eq.0) write(*,*)"Warning for",k," filter: Nused=0"
      if(nused.ne.0)then
         f=chi/dble(nused)
      else
         f=1.d9
      endif

      end




