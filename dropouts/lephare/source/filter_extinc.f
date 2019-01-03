c
      program ext_curv
c
c     Compute the mean extinction for a given set of filter 
c
      implicit none
      integer*4 nbf,maxsize,wmax
      parameter (maxsize=110000)
c      parameter (wmax=8000)
      INCLUDE 'dim_wave.decl'      
      INCLUDE 'dim_filt.decl'
      integer*4 test,lnblnk,iext,ifilt,i,j,imag,iextg
      integer*4 jmax(nbf),id(nbf),calib(nbf)
      real*8 lext(wmax),ext(wmax)
      real*8 lextg(wmax),extg(wmax)
      real*8 lfilt(wmax),filt(wmax),aint,albd,albdav,rv
      character*4096  zpdir,zpwork,paravc(500)
      character*4096  file,filters,extc,output,galc
      character*4096  name(nbf),param
      character*4096  str
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c  environment variable 
      call getenv('LEPHAREDIR',zpdir)
      test=lnblnk(zpdir)
      if (test .eq. 0) then
        write(6,*) 'WARNING :  variable LEPHAREDIR not defined'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(6,*) 'WARNING :  variable LEPHAREWORK not defined'
        stop
      endif
cccccccccccccccccccccccccc
c help option
      param='filter_extinc'
      call get_help(test)
      if (test .eq. 1) call help(param)      
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c Initialisation of Input parameter from config file
c  read option
      param='-FILTER_FILE'
      call opt_line(param,1,paravc,test)
      if (test.eq.1)  filters=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  then
          param='-f'
          call opt_line(param,1,paravc,test)
      endif
      if (test.eq.1)  filters=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1) call err_option(param,1)
      param='-EXT_CURV'
      call opt_line(param,1,paravc,test)
      if (test.eq.1)  extc=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  then
          param='-e'
          call opt_line(param,1,paravc,test)
      endif
      if (test.eq.1)  extc=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  then 
         call err_option(param,2) 
         extc='NONE'
      endif
      param='-GAL_CURV'
      call opt_line(param,1,paravc,test)
      if (test.eq.1)  galc=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  then
          param='-g'
          call opt_line(param,1,paravc,test)
      endif
      if (test.eq.1)  galc=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  then 
         call err_option(param,2) 
         galc='CARDELLI'
      endif
      param='-OUTPUT'
      call opt_line(param,1,paravc,test)
      if (test.eq.1)  output=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  then
          param='-o'
          call opt_line(param,1,paravc,test)          
      endif
      if (test.eq.1)  output=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  then 
        call err_option(param,2)
        output='NONE' 
      endif
      if (output .ne. 'NONE') open(2,file=output,status='unknown') 
cccccccccccccccccccccccccccc
c  writing options 
       if (output .ne. 'NONE') 
     > write(2,'(A)')"# Computing ATMOSPHERIC AND GALACTIC EXTINCTION #"
       write(6,'(A)')"########################################"
       write(6,'(A)')"# Computing the ATMOSPHERIC/GALACTIC EXTINCTION #"
       write(6,'(A)')"#     with the following OPTIONS       #"
       file = zpwork(1:lnblnk(zpwork)) //'/filt/'
     >      // filters(1:lnblnk(filters))
       write(6,'(2A)')"# FILTER_FILE : ",file(1:lnblnk(file))
       if (output .ne. 'NONE')
     > write(2,'(2A)')"# FILTER_FILE : ",file(1:lnblnk(file))
       if (extc .ne. 'NONE') then 
          file = zpdir(1:lnblnk(zpdir)) //'/ext/'
     >      // extc(1:lnblnk(extc))
       else
          file = 'NONE'
       endif
       write(6,'(2A)')"# EXT_CURVE   : ",file(1:lnblnk(file))  
       if (output .ne. 'NONE')
     > write(2,'(2A)')"# EXT_CURVE   : ",file(1:lnblnk(file))

       if (galc .ne. 'CARDELLI') then 
          file = zpdir(1:lnblnk(zpdir)) //'/ext/'
     >      // galc(1:lnblnk(galc))
       else
          file = 'CARDELLI law'
       endif
       write(6,'(2A)')"# GAL_CURVE   : ",file(1:lnblnk(file))  
       if (output .ne. 'NONE')
     > write(2,'(2A)')"# GAL_CURVE   : ",file(1:lnblnk(file))  

       write(6,'(2A)')"# OUTPUT      : ",output(1:lnblnk(output))  
       if (output .ne. 'NONE')
     > write(2,'(2A)')"# OUTPUT      : ",output(1:lnblnk(output))  
       write(6,'(A)')"###########################################"
       write(6,'(A)')" Filters Ext(mag/airmass) Albda/Av Albda/E(B-V) " 
       if (output .ne. 'NONE')
     > write(2,'(A)')"#Filters Ext(mag/airmass) Albda/Av Albda/E(B-V) " 
c  read atmospheric extinction 
      if (extc .ne. 'NONE') then
         file = zpdir(1:lnblnk(zpdir)) //'/ext/'
     >      // extc(1:lnblnk(extc))  
         open(1,file=file,status='unknown')
         i=0
         do while (.true.)
            i=i+1
	    read(1,*,end=10)  lext(i), ext(i)
         enddo
 10      iext=i-1  
         close(1)
      endif
c  read galactic extinction 
      if (galc .ne. 'CARDELLI') then
         file = zpdir(1:lnblnk(zpdir)) //'/ext/'
     >      // galc(1:lnblnk(galc))  
         open(1,file=file,status='unknown')
         i=0
         do while (.true.)
            i=i+1
	    read(1,*,end=11)  lextg(i), extg(i)
         enddo
 11      iextg=i-1  
         close(1)
      endif
c  read the filters from file and write results in output
c
      file = zpwork(1:lnblnk(zpwork)) //'/filt/'
     >        // filters(1:lnblnk(filters))
      open(1,file=file,status='unknown',err=56)
      read(1,'(A)') str 
      call val_string(str,paravc,test)
      read(paravc(2),'(i10)') imag
c  check Rv to be applied 
      rv = 3.1
      if (galc .ne. 'CARDELLI') then
         if ( galc .eq. 'SMC_prevot.dat') then
            rv=2.72
         elseif  ( galc(1:7) .eq. 'SB_calz'
     >        .or. galc(1:6) .eq. 'calzet') then
            rv=4.05
         else
            rv=3.1
         endif
      else
         rv=3.1
      endif
      if (output .eq. 'NONE')  then
        write(6,'(A,2x,f6.2,2x,A)') 
     >             'assuming Rv=',rv,' for this Extinction law '
      endif
c
      do i =  1, imag
         read(1,'(A)') str
         call val_string(str,paravc,test)
         read(paravc(2),'(i10)') jmax(i)
           read(paravc(4),'(i10)') calib(i)
           read(paravc(5),'(i10)') id(i)
         name(i)=paravc(3)
         do j = 1,jmax(i)
	    read(1,*)  lfilt(j), filt(j) 
 	 enddo
         ifilt=jmax(i)
c  compute atmospheric extinction 
       if (extc .ne. 'NONE') then
         call comp_ext(lfilt,filt,ifilt,lext,ext,iext,aint)
       else
         aint = 99.99
       endif
c  compute galactic extinction 
       if (galc .ne. 'CARDELLI') then
         call comp_ext(lfilt,filt,ifilt,lextg,extg,iextg,albd)
c galactic curves given in k(lbda) (=A(lbda)/E(B-V))
c -> A(lbda)/Av = A(lbda)/E(B-V) / Rv)
c Rv=3.1 except for Calzetti law (4.05) and SMC Prevot (2.72)
         if ( galc .eq. 'SMC_prevot.dat') then
            albdav = albd/2.72
            rv=2.72
         elseif  ( galc(1:7) .eq. 'SB_calz'
     >        .or. galc(1:6) .eq. 'calzet') then
            albdav = albd/4.05
            rv=4.05
         else
            albdav = albd/3.1
            rv=3.1
         endif
       else
         call comp_gal(lfilt,filt,ifilt,albdav)
c  output A(lbd)/Av->A(lbd)/(Rv*E(B-V))->A(lbd)/E(B-V)=Rv*A(lbd)/Av
         albd = albdav*3.1
         rv=3.1
       endif
       if (output .ne. 'NONE')          
     >    write(2,110) name(i),aint,albdav,albd
         write(6,110) name(i),aint,albdav,albd
      enddo
      close(1) 
      if (output .ne. 'NONE') close(2)
c
 110  format(A10,2x,3(f8.3,2x))
      stop
 56   write (6,*) 'File ',file(1:lnblnk(file)),
     >  ' not found -> STOP '
   
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine comp_ext(lfilt,filt,ifilt,lext,ext,iext,aint)
      implicit none 
      integer*4 wmax,ifilt,i,iext,smax
c      parameter (wmax=8000)
      INCLUDE 'dim_wave.decl'      
      real*8 lext(wmax),ext(wmax)
      real*8 lfilt(wmax),filt(wmax)
      real*8 ls(wmax),fs(wmax),ts(wmax),aint,fint
c
c  sampling the curves
      call sampling(lext,ext,iext,lfilt,filt,ifilt,ls,fs,ts,smax)   
c
c   integrate the extinction curve through the filter
      fint = 0
      aint = 0 
      do i = 1, smax -1
         fint= fint+( ts(i) + ts(i+1) ) /2.  * (ls(i+1) - ls(i))
         aint= aint+(ts(i)+ts(i+1))*(ls(i+1)-ls(i))*(fs(i)+fs(i+1))/4.
      enddo
      aint = aint / fint 
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine comp_gal(lfilt,filt,ifilt,acor)
      implicit none 
      integer*4 wmax,ifilt,i,iextg,smax
c      parameter (wmax=8000)
      INCLUDE 'dim_wave.decl'      
      real*8 lextg(wmax),extg(wmax)
      real*8 lfilt(wmax),filt(wmax),lmin,lmax
      real*8 ls(wmax),fs(wmax),ts(wmax),acor,fint
c
      lmin = lfilt(1)
      lmax = lfilt(ifilt)
c
c  computes the galactic extinction 
      call gal_curv(lmin,lmax,lextg,extg,iextg)
c
c  sampling the curves
      call sampling(lextg,extg,iextg,lfilt,filt,ifilt,ls,fs,ts,smax)   
c
c   integrate the extinction curve through the filter
      fint = 0
      acor = 0 
      do i = 1, smax -1
         fint =fint+( ts(i) + ts(i+1) ) /2.  * (ls(i+1) - ls(i))
         acor =acor+(ts(i)+ts(i+1))*(ls(i+1)-ls(i))*(fs(i)+fs(i+1))/4.
      enddo
      acor = acor / fint 
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Computes the galactic extinction based on Cardelli et al., 1989, ApJ 345
      subroutine gal_curv(lmin,lmax,lextg,extg,iextg)
      implicit none 
      integer*4 wmax,i,iextg
       INCLUDE 'dim_wave.decl'      
c     parameter (wmax=8000)
      real*8 lmin,lmax,lextg(wmax),extg(wmax)
      real*8  dlbd,ext_gal
      external ext_gal
c
c  step in lbda 
      dlbd = (lmax-lmin)/400. 
c
      do i = 1, 402
         lextg(i) = lmin + (i-1)*dlbd
         extg(i) = ext_gal(lextg(i))
      enddo   
      iextg = 402
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccc
c  compute albd/av at a given lambda (A) 
      REAL*8 FUNCTION ext_gal(lb)
c
      implicit none 
c
      real*8 lb,rv,f1,f2,x,y,fa,fb
c
      rv=3.1
      x = 10000./lb
      y = x - 1.82
      f1=0
      f2=0
      if (x .le. 1.1) then
         f1 =  0.574 * x**(1.61)
         f2 = -0.527 * x**(1.61)
      elseif ( x.gt. 1.1 .and. x.lt. 3.3 ) then 
          f1 = 1+0.17699*y -0.50447*y**2 -0.02427*y**3 +0.72085*y**4
          f1 = f1+0.01979*y**5 -0.77530*y**6 +0.32999*y**7 
          f2 = 1.41338*y +2.28305*y**2 +1.07233*y**3
          f2 = f2-5.38434*y**4-0.62251*y**5+5.30260*y**6 -2.09002*y**7
      elseif ( x.ge. 3.3 .and. x.lt. 5.9 ) then 
          f1 = 1.752 -0.316*x -0.104/((x-0.467)**2+0.341) 
          f2 = -3.090 +1.825*x +1.206/((x-4.62)**2+0.262)  
      elseif ( x.ge. 5.9 .and. x.lt. 8 ) then 
          fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
          fb = 0.2130*(x-5.9)**2 + 0.1207*(x-5.9)**3 
          f1 =  1.752 -0.316*x -0.104/((x-0.467)**2+0.341) +fa 
          f2 = -3.090 +1.825*x +1.206/((x-4.62)**2+0.262)  +fb
      elseif ( x.ge. 8 ) then 
          f1 = -1.073 -0.628*(x-8)  +0.137*(x-8)**2 -0.070*(x-8)**3 
          f2 =  13.670 +4.257*(x-8)  -0.420*(x-8)**2 +0.374*(x-8)**3 
      endif
      ext_gal = f1 + f2/rv 
c
      return
      end
c

