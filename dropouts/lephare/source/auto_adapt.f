ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     21/07/04
c     subroutine to minimize the differences between the predicted 
c     apparent magnitudes and the observed.
c     The minimization is done using the minuit library.
      subroutine auto_adapt(chiin,chitot,a0in,a0fit,
     .      a1in,a1fit,a2in,a2fit,a3in,a3fit,realise,degre,adcont,
     .      residu,adapterror)

      implicit none

      integer*4  i,k,imagm,nbf,realise,fl1,fl2,fl_auto
      integer*4  meth_ada,zadapt,bdincl
      INCLUDE 'dim_filt.decl'
      parameter (zadapt=30000)
      character*4096 adapterror,zp_med,adaptMed
      real*8     cont_ada(zadapt),adcont
      integer*4  mod_ada(zadapt)
      integer*4  degre,nused(imagm),nb,ngals_ada,bused(nbf)
      real*8 ab_ada(zadapt,nbf),sab_ada(zadapt,nbf)
      real*8 magm_ada(zadapt,nbf),zs_ada(zadapt)
      real*8 a0fit(imagm),a0in(imagm)
      real*8 a1fit(imagm),a1in(imagm)
      real*8 a2fit(imagm),a2in(imagm)
      real*8 a3fit(imagm),a3in(imagm)
      real*8 residu,residuc,chiin,chifit,chitot,sig,sig2(imagm)
      real*8 min_err(500),diffminerr,tot,av_error,fac_err
      integer*4 ind_med(zadapt)
      real*8 medDiff(imagm),med_mag(zadapt)
      real*8 med_mag2(zadapt)

      real*8 vstrt(4),stp(4),zero,u(500),gin
      integer*4 nprm(4),isw(7),ierflg
      character*4096 pnam(4)
      data nprm /   1   ,    2 , 3 , 4  /
      data pnam /'a0'   ,   'a1'  , 'a2' , 'a3'/
      data vstrt/   0.  ,    0.  , 0., 0./
      data stp  /   0.02 ,    0.02 ,  0.02,  0.02/
      data zero / 0./      

      external fcnfunc
      external bdincl 

c     common for minuit
      common/func_int/k,ngals_ada,fl1,fl2,fl_auto,meth_ada,imagm
      common/func_real/sig
      common /func_tab/ ab_ada,sab_ada,magm_ada,zs_ada
      common /func_tab2/ cont_ada,mod_ada
      common /min_err/min_err
      common /fac_err/fac_err
      common /zp_med/adaptMed

      common/mn7ext/u
      common/mn7flg/isw

      chitot=0.d0


c     Loop on each filter
      do k=1,imagm

c      Initialization       
       a0fit(k)=0.
       a1fit(k)=0.
       a2fit(k)=0.
       a3fit(k)=0.
    
       do i=1,zadapt
         med_mag(i)=0.d0
       enddo
       medDiff(k)=0

       if(bdincl(k-1,adcont,imagm-1).eq.1.or.adcont.le.0)then

c       Compute sigma for 3 sigma clipping and check if enougth gal
c       Compute the median of Obs - predicted
        nused(k)=0
        sig=0.d0
c       Loop over objets choose for adapt
        do i=1,ngals_ada 
c         Check the filter
          bused(k)=bdincl(k-1,cont_ada(i),imagm-1)
          if(cont_ada(i).eq.0)bused(k)=1
          if(bused(k).eq.1)then
c          If the mag can be used
           if(dabs(ab_ada(i,k)).le.80.d0.and.sab_ada(i,k).gt.0.d0)then
            if(dabs(magm_ada(i,k)).le.80.d0)then
c               Compute the sigma for the 3 sigma clipping
                sig=sig+(ab_ada(i,k)-magm_ada(i,k))**2.
                nused(k)=nused(k)+1
c               Compute the diff for the median
c               I remove the previous correction from the difference
                med_mag(nused(k))=-magm_ada(i,k)+a0in(k)+ab_ada(i,k)
c               Compute the |diff| for the 1.48*median (~sigma)
c               I remove the previous correction from the difference
                med_mag2(nused(k))=abs
     .                            (-magm_ada(i,k)+a0in(k)+ab_ada(i,k))
            endif
           endif
          endif
        enddo
        if(nused(k).gt.1)then
          sig=dsqrt(sig/(nused(k)-1))
          sig2(k)=sig
          call sortNR(med_mag,nused(k),zadapt,ind_med)
          medDiff(k)=med_mag(ind_med(nint(dble(nused(k))/2.)))
          call sortNR(med_mag2,nused(k),zadapt,ind_med)
          sig2(k)=med_mag2(ind_med(nint(dble(nused(k))/2.)))       
        else
          sig=999.9
        endif    
c        write(*,*)"sigma for ",k,", :",sig," and median ",
c     .   medDiff(k),"with ",nused(k)," obj"
cccccccccccccccccccccccc
        write(41,*)"Filter :",k
c       check if the number of objects allow the minimization 
        if(nused(k).ne.0)then
c         write minuit output in minuit.dat
          call mninit(5,41,7)
c         initialisation of para
          do i= 1,4
           call mnparm(nprm(i),pnam(i),vstrt(i),stp(i),zero,zero,ierflg)
           if (ierflg .ne. 0)  then
            write (6,'(a,i3)')  ' unable to define parameter no.',i
            stop
           endif
          enddo
          call mnseti('Minimization mag auto - library')
          call mncomd(fcnfunc,'set str 2',ierflg,0)
c         If you want to fix one free parameter (degre of the poly for the fit)
          if(degre.lt.1)call mncomd(fcnfunc,'fix 1',ierflg,0)
          if(degre.lt.2)call mncomd(fcnfunc,'fix 2',ierflg,0)
          if(degre.lt.3)call mncomd(fcnfunc,'fix 3',ierflg,0)
          if(degre.lt.4)call mncomd(fcnfunc,'fix 4',ierflg,0)
c         put limit for the free parameters
          call mncomd(fcnfunc,'set lim 1 -1. 1.',ierflg,0)
          call mncomd(fcnfunc,'set lim 2 -1. 1.',ierflg,0)
          call mncomd(fcnfunc,'set lim 3 -1. 1.',ierflg,0)
          call mncomd(fcnfunc,'set lim 4 -1. 1.',ierflg,0)
c         derived hessian matrix
          call mnexcm(fcnfunc,'hes',zero,0,ierflg,0)
c         migrad to converge
          call mnexcm(fcnfunc,'migrad',zero,0,ierflg,0)
c         errors
          call mnexcm(fcnfunc,'minos',zero,0,ierflg,0)
c         write the solution
          call mnexcm(fcnfunc,'sho par',zero,0,ierflg,0)
c         chi2 of the fit
          call fcnfunc(4,gin,chifit,u,4)
c         save values of the fit
          a0fit(k)=u(1)
          a1fit(k)=u(2)
          a2fit(k)=u(3)
          a3fit(k)=u(4)
          chitot=chitot+chifit

        else

c        no adaptation 
         a0fit(k)=0.d0
         a1fit(k)=0.d0
         a2fit(k)=0.d0
         a3fit(k)=0.d0

        endif
       
       endif

      enddo

c     New coeficients of the poly function is the sum of old and last fit
      residuc=0.0d0
      tot=0.d0
      do k=1,imagm
c      CORRECTION !
       if(a0fit(k).gt.0.d0) residuc=residuc+abs(a0fit(k))
ccccccccccccccccccc
       a0fit(k)=a0fit(k)+a0in(k)
       a1fit(k)=a1fit(k)+a1in(k)
       a2fit(k)=a2fit(k)+a2in(k)
       a3fit(k)=a3fit(k)+a3in(k)
       tot=tot+a0fit(k)
      enddo
ccccccccccccccccccccc
c      Average of the coef at 0
      do k=1,imagm
cc       a0fit(k)=a0fit(k)-tot/dble(imagm)
c       if(k.ne.fl_auto)a0fit(k)=a0fit(k)-a0fit(fl_auto)
c       medDiff(k)=medDiff(k)-medDiff(fl_auto)
       if (adaptMed(1:1).eq.'Y'.or.adaptMed(1:1).eq.'y')then 
         a0fit(k)=medDiff(k) 
c         write(25,*)k,a0fit(k),medDiff(k),a0in(k)
       endif
      enddo      
c      a0fit(fl_auto)=0
ccccccccccccccccccccc


cc    condition to stop the iteration
      if(chitot.eq.0)then
       residu=0
      else
       residu=abs(chitot-chiin)/chitot
      endif
c      write(6,*) residu
      write(6,'(A,1x,I6,1x,A,E12.6,1x,A)') "  ",ngals_ada,"objects:
     > chi2 variation over all filters of ",residu*100,"%"
      write(6,'(A,1x,I6,1x,A,f8.4,1x,A)') "  ",ngals_ada,"objects:
     > average zero-point variation of ",residuc/dble(imagm)," mag"
c     Stop if all the variation are less than 0.03 mag
      realise=2
      do k=1,imagm
        if(abs(a0fit(k)-a0in(k)).gt.0.03)then
          realise=0
        endif
      enddo
c     Stop if the residual chi2 varues less than 5%
c      if(residuc/dble(imagm).le.1.d-3)realise=2
      if(residu.le.1.d-2)realise=2

c      Tentative pour adapter les erreurs
c      Compute the error scale to apply
c      Loop over the filter
       do k=1,imagm
        nused(k)=0
        sig=0.d0
        av_error=0.
        nb=0
        if(bdincl(k-1,adcont,imagm-1).eq.1.or.adcont.le.0)then
c        Loop over objets choose for adapt
         do i=1,ngals_ada 
          bused(k)=bdincl(k-1,cont_ada(i),imagm-1)
          if(bused(k).eq.1)then
           if(dabs(ab_ada(i,k)).le.80.d0.and.sab_ada(i,k).gt.0.d0)then
            if(dabs(magm_ada(i,k)).le.80.d0)then
c            Observed - ( predicted in - olf zp + new zp) 
             diffminerr=ab_ada(i,k)-magm_ada(i,k)+a0in(k)-a0fit(k)
c     .                   +a0fit(fl_auto)
c             write(25,*)k,i,diffminerr,ab_ada(i,k)-magm_ada(i,k)
c             if(diffminerr.le.3.*sig2(k))then
              nused(k)=nused(k)+1
              nb=nb+1
              sig=sig+diffminerr**2.
              av_error=sab_ada(i,k)+av_error
c             Compute the |diff| for the 1.48*median (~sigma)
              med_mag2(nused(k))=abs(diffminerr)
c             median of the error
              med_mag(nused(k))=sab_ada(i,k)
c             endif
            endif
           endif
          endif
         enddo
         if(nused(k).gt.1)then
           sig=dsqrt(sig/dble(nb))
           av_error=av_error/dble(nb)
           call sortNR(med_mag2,nused(k),zadapt,ind_med)
c           write(*,*)"sig",k,sig,
c     .               med_mag2(ind_med(nint(dble(nused(k))/2.)))       
           sig=med_mag2(ind_med(nint(dble(nused(k))/2.)))
           call sortNR(med_mag,nused(k),zadapt,ind_med)
c           write(*,*)"med",k,av_error,
c     .                med_mag(ind_med(nint(dble(nused(k))/2.)))       
           av_error=med_mag(ind_med(nint(dble(nused(k))/2.))) 
           if((adapterror(1:1).eq."Y".or.adapterror(1:1).eq."y"))then 
            if(av_error.lt.sig)then
             min_err(k)=dsqrt(  sig**2. - av_error**2.)/fac_err
             min_err(k)=max(min_err(k),0.01)
            else
             min_err(k)=0.01
            endif
           endif    
         endif
        endif
      write
     .(*,'(a,(f8.6,1x),a,1(f11.6,1x),a,1(i5,1x),3(a,1(f11.6,1x)))')
     . "auto-adapt :",a0fit(k)," sigma : ",sig, " nb obj : ",nused(k),
     ."av error ",av_error, "min error",min_err(k)," Median ",medDiff(k)
       enddo
c      endif
c      endif


c10000 format("a0,a1,object ->",2(f8.5,1x),I7,a1,$)

      end

