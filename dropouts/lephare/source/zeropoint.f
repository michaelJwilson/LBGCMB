c      last modif 06/01/2000
c      Author : S. ARNOUTS 
c      description : compute the Vega Zeropoint and AB correction
c       for different filters 
c      and gives the characteristics of filters 
c      Input : -Vega spectrum from Fioc et al., 1999 
c              -Filter list : filter_xxx.dat: Arnouts et al., 1999
c      subroutines USED :  sampling.f  
c      Environnement variable : LEPHAREDIR local directory of LEPHARE
c      to be defined with setenv LEPHAREDIR local_dir_zphot 
      subroutine zeropoint(filtfile,veg,ab,lmoy,width,fcorr,imag)
      implicit none
      integer*4 nbf,maxsize,wmax
      INCLUDE 'out_unit.decl' 
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_wave.decl'
      parameter (maxsize=110000)
      integer*4 i,imag,j,k,jmax(nbf),id(nbf),calib(nbf)
      integer*4 test,pass,smax,imax1,imax2,imax0,imaxbd,imaxsun
      integer*4 lnblnk
      real*8    c,lmean,trans,lmoy(nbf),area(nbf),arean(nbf)
      real*8    wav(5000), fl(5000),fpic(nbf),lbpic(nbf)
      real*8    wavbd(5000), flbd(5000), fmeltg(nbf)
      real*8    lbsup(nbf),lbinf(nbf),width(nbf),frac,slope
c      real*8    lamb(nbf,maxsize),rep(nbf,maxsize),conv,tmax
      real*8    lamb(nbf,wmax),rep(nbf,wmax),conv,tmax
      real*8    fmenu(nbf), fmel(nbf), ab(nbf),veg(nbf),tg(nbf)
      real*8    lnum(nbf),lden(nbf),leff(nbf),fvega(nbf)
c      real*8    ls(maxsize),fs(maxsize),ts(maxsize)
c      real*8    x2(maxsize),y2(maxsize),flux,dlbd
      real*8    ls(wmax),fs(wmax),ts(wmax)
      real*8    x2(wmax),y2(wmax),flux,dlbd
      real*8    wavsun(5000),lsun(5000),msun(nbf)
      real*8    num0,den0,num1,den1,num2,den2,num3
      real*8    leff0(nbf),leff1(nbf),leff2(nbf)
c      real*8    leff3(nbf)
c      real*8    corr0(nbf),corr1(nbf),corr2(nbf),corr3(nbf)
      real*8    hckt,bb,x0(1001),y0(1001),dlbd0,bb2
      real*8    l_eff(nbf),fcorr(nbf)
      parameter(c=2.99792458e+18)
      character*4096  name(nbf)
      character*4096  vega,sun,filtfile,filtfile2,bd17
      character*4096  zpdir,zpwork,paravc(500)
      character*4096 str

c     If you want to give a name to the output screen file
      if(UO.eq.30)
     .      open(UO,file='/tmp/screenZeroPoint.dat',status='unknown') 
c
c
c  environmental variable 
      call getenv('LEPHAREDIR',zpdir)
      test=lnblnk(zpdir)
      if (test .eq. 0) then
        write(UO,*) 'WARNING :  variable LEPHAREDIR not defined'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(UO,*) 'WARNING :  variable LEPHAREWORK not defined'
        stop
      endif
c  input files  
      vega = zpdir(1:lnblnk(zpdir)) // '/vega/VegaLCB.sed'  
      bd17 = zpdir(1:lnblnk(zpdir)) // '/vega/BD+17o4708.sed'
      sun  = zpdir(1:lnblnk(zpdir)) // '/vega/SunLCB.sed'
c
      filtfile2 = zpwork(1:lnblnk(zpwork)) // '/filt/' 
     >         //filtfile(1:lnblnk(filtfile)) 
cccccccccccccccccccccccccccccccccccccccccccccccc 
c  reading the Vega spectrum 
c  use vegaLCB.sed file for Vega spectrum from Fioc et al., 1999
        open(1, file=vega,status='old')
        read(1,*)   
        i=0
	do while (.true.)
	  i=i+1
	  read(1,*,end=10)  wav(i), fl(i)
 	enddo
 10     imax1=i-1  
        close(1)
cccccccccccccccccccccccccccccccccccccccccccccccc
c  use BD+17o4708.sed file for Thuan Gunn 
        open(1, file=bd17,status='old')
        read(1,*)   
        i=0
	do while (.true.)
	  i=i+1
	  read(1,*,end=11)  wavbd(i), flbd(i)
 	enddo
 11     imaxbd=i-1  
        close(1)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c  use Sun Luminosity from Fioc et al., 1999
        open(1, file=sun,status='old')
        read(1,*)   
        i=0
	do while (.true.)
	  i=i+1
	  read(1,*,end=12)  wavsun(i), lsun(i)
 	enddo
 12     imaxsun=i-1  
        close(1) 
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c  reading name, lambda, T(lbda) of filters and defines the WIDTH
        open(1,file=filtfile2,status='old',err=56)
        read(1,'(A)') str
        call val_string(str,paravc,test)
        read(paravc(2),'(i10)') imag
        do i =  1, imag
           read(1,'(A)') str
           call val_string(str,paravc,test)
           read(paravc(2),'(i10)') jmax(i)
           name(i) = paravc(3)
           read(paravc(4),'(i10)') calib(i)
           read(paravc(5),'(i10)') id(i)
           fpic(i)=0
           do j = 1,jmax(i)
	      read(1,*)  lamb(i,j), rep(i,j)
c  look at the peak of the filter 
              if (rep(i,j).ge.fpic(i) ) then
                  fpic(i) =  rep(i,j)
                  lbpic(i)= lamb(i,j)
              endif
 	   enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  define the width based on a fraction of the Peak level
           frac=0.5
           pass=0
           do j = 1,jmax(i)      
              if (rep(i,j).le.(fpic(i)*frac) .and. 
     >            rep(i,j+1).gt.(fpic(i)*frac) .and. pass.eq.0 ) then
c  No resampling if Delta(lambda)<5A
                  if ((lamb(i,j+1)-lamb(i,j)).le.5) then 
                     lbinf(i)=  (lamb(i,j+1)+lamb(i,j))/2.
                     pass=1
                  else
c Resampling Delta(lambda) with 1A
                     slope=(rep(i,j+1)-rep(i,j))/
     >                      (lamb(i,j+1)-lamb(i,j))
                     do k=1,IDNINT(lamb(i,j+1)-lamb(i,j))
                       if ((rep(i,j)+slope*float(k)).ge.
     >                      (fpic(i)*frac).and. pass.eq.0 ) then
                          lbinf(i)= lamb(i,j)+DBLE(k) 
                          pass=1
                       endif
                     enddo
                  endif
               elseif   (rep(i,j).ge.(fpic(i)*frac) .and.
     >            rep(i,j+1).lt.(fpic(i)*frac)) then
c  No resampling if Delta(lambda)<5A
                  if ((lamb(i,j+1)-lamb(i,j)).le.5) then 
                     lbsup(i)=  (lamb(i,j+1)+lamb(i,j))/2.
                  else
c Resampling Delta(lambda) with 1A
                     slope=(rep(i,j+1)-rep(i,j))/
     >                      (lamb(i,j+1)-lamb(i,j))
                     do k=1,IDNINT(lamb(i,j+1)-lamb(i,j))
                       if ((rep(i,j)+slope*DBLE(k)).le.
     >                      (fpic(i)*frac) ) then
                          lbsup(i)= lamb(i,j)+float(k) 
                         goto 63
                       endif
                     enddo
 63                  continue
                  endif
               endif
           enddo   ! j
           width(i)=lbsup(i)-lbinf(i)      
        enddo
        close(1) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  computation AB mag , characteristics  for each filter
      write(6,'(A)') '# NAME    IDENT      Lbda_mean    Lbeff(Vega)  
     >  FWHM     AB-cor  TG-cor  VEGA  M_sun(AB) CALIB    Lb_eff   Fac_c
     >orr'
      do i = 1, imag
        fmel(i)  = 0.
        fmeltg(i)= 0.
        ab(i)    = 0.
        veg(i)   = 0.  
        tg(i)    = 0.
        fmenu(i) = 0.
        area(i)  = 0.
        arean(i) = 0.
        lnum(i)  = 0.  
        lden(i)  = 0.  
        lmoy(i)  = 0.
        fvega(i)=0
c   REsampling the 2 curves with same absisse
        do j = 1, jmax(i) 
           x2(j) = lamb(i,j)
           y2(j) = rep(i,j)
        enddo
        imax2 = jmax(i)
        call sampling(wav,fl,imax1,x2,y2,imax2,ls,fs,ts,smax)  
        tmax=0
        do j = 1,smax 
           if (ts(j).ge.tmax) tmax=ts(j)
        enddo   
        if (x2(1).ge.wav(1).and.x2(jmax(i)).le.wav(imax1))then
          do j = 1,smax - 1
            lmean   = (ls(j) + ls(j+1))/2.
            trans   = (ts(j) + ts(j+1))/(2.*tmax)
            flux    = (fs(j)+fs(j+1))/2.
            dlbd    =  ls(j+1)-ls(j) 
            conv    = c/lmean**2        
            lmoy(i) = lmoy(i) + lmean*trans*dlbd
            area(i) = area(i) + trans*dlbd
            lnum(i) = lnum(i)+flux*trans*lmean*dlbd
            lden(i) = lden(i)+flux*trans*dlbd
            arean(i)= arean(i) + trans*dlbd*conv
            fmel(i) = fmel(i) + flux*trans*dlbd  
            fvega(i)= fvega(i)+ flux*dlbd
          enddo
          lmoy(i) = lmoy(i) / area(i)
          ab(i)   = -2.5*dlog10(fmel(i)/arean(i)) -48.59
          veg(i)  = +2.5*dlog10(fmel(i)/area(i))  
          fvega(i)= -2.5*dlog10(fvega(i))
          leff(i) = lnum(i) / lden(i)
        else
          do j = 1,smax - 1
            lmean   = (ls(j) + ls(j+1))/2.
            trans   = (ts(j) + ts(j+1))/(2.*tmax)
            dlbd    =  ls(j+1)-ls(j) 
            area(i) = area(i) + trans*dlbd
            lmoy(i) = lmoy(i) + lmean*trans*dlbd
          enddo
          lmoy(i) =  lmoy(i) / area(i)
          ab(i)   = -99.99
          veg(i)  = -99.99
          fvega(i)= -99.99
          leff(i) = -99.99*10000 
        endif 
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c  other calibration for Lambda eff 
        do j = 1, jmax(i) 
           x2(j) = lamb(i,j)
           y2(j) = rep(i,j)
        enddo
        imax2 = jmax(i)
c   build a 1000 elements sampling inside each filter
        dlbd0 = (x2(jmax(i))-x2(1))/1000.
        imax0 = 1001
        do j = 1, 1001
          x0(j) = x2(1) + (j-1)*dlbd0
          y0(j) = 0
        enddo  
c       write(*,*) x2(1),x2(jmax(i)),dlbd0,x0(1),x0(1001) 
        call sampling(x0,y0,imax0,x2,y2,imax2,ls,fs,ts,smax)   
c   constant hc/KT for Black Body with T=10000K
c    Blbda= 2.h.c^2/lbda^5/(exp(hc/kTlbda)-1)
        num0=0.
        den0=0.
        num1=0.
        den1=0.
        num2=0.
        den2=0.
        num3=0.
c               
        hckt    = 14397.64648         
        do j = 1,smax-1
          lmean   = (ls(j) + ls(j+1))/2.
          trans   = (ts(j) + ts(j+1))/2.
          dlbd    =  ls(j+1)-ls(j) 
          conv    = c/lmean**2    
          if ( (dexp(hckt/lmean)-1) .lt. 1.e-5) then 
             bb = 1./(hckt/lmean)/lmean**5
          else
            bb= 1./(dexp(hckt/lmean)-1)/lmean**5
          endif
          num3    = num3 + lmean*trans*dlbd
c    nuBnu=Cte -> Blbda=1/lbda
          num0    = num0 + trans*dlbd
          den0    = den0 + trans*dlbd/lmean
c    Bnu=nu -> Blbda=1/lbda^3
          num1    = num1 + trans/lmean/lmean*dlbd
          den1    = den1 + trans/lmean/lmean/lmean*dlbd
c    Blbda= BB(10K)
          num2    = num2 + trans*bb*lmean*dlbd
          den2    = den2 + trans*bb*dlbd
        enddo
c  Lmean  for Bnu=ctt
        lmoy(i) = num3/num0 
c  Leff for nuBnu=cte
        leff0(i) = num0/den0
c  Leff for Bnu=nu
        leff1(i) = num1/den1
c  Leff for BB
        leff2(i) = num2/den2        
c  compiles values according to FILTER_CALIB  
c  <F>lp,cor = <F>lp * int(Tn dn) /  int(Bn/Bno Tn dn)     
c     fcorr' = int(Bn/Bno Tn dn) / int(Tn dn)
c  <F>lp,cor = <F>lp * flcorr 
c  with  fcorr=1/fcorr'. It is measured at the end. 
        if (calib(i).le.0 .or. calib(i).gt.6) then
           l_eff(i) =  lmoy(i)
           fcorr(i) = 1.0
        elseif  (calib(i).eq.1) then
           l_eff(i) =  leff0(i)
           fcorr(i) =  1/leff0(i)*den0/num1
        elseif  (calib(i).eq.2) then
           l_eff(i) =  leff1(i)
           fcorr(i) =  leff1(i)*den1/num1
        elseif  (calib(i).eq.3) then
           l_eff(i) =  leff2(i)
c           computing BB at lbd_eff derived with BB
           if ( (dexp(hckt/leff2(i))-1) .lt. 1.e-5) then 
              bb2 = 1./(hckt/leff2(i))/leff2(i)**5
           else
              bb2 = 1./(dexp(hckt/leff2(i))-1)/leff2(i)**5
           endif
           fcorr(i) =  1./leff2(i)/leff2(i)/bb2*den2/num1 
        elseif  (calib(i).eq.4) then
           l_eff(i) =  leff0(i)
c           computing BB at lbd_eff derived with nuBnu=ctt
           if ( (dexp(hckt/leff0(i))-1) .lt. 1.e-5) then 
              bb2 = 1./(hckt/leff0(i))/leff0(i)**5
           else
              bb2 = 1./(dexp(hckt/leff0(i))-1)/leff0(i)**5
           endif 
           fcorr(i) = 1./leff0(i)/leff0(i)/bb2*den2/num1 
        elseif  (calib(i).eq.5) then
           l_eff(i) =  leff0(i)
           fcorr(i) =  leff0(i)*den1/num1
        endif
c  we apply the reverse as correction : fcorr =  int(Tn dn)/ int(Bn/Bno Tn dnu)
         fcorr(i) = 1/fcorr(i)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  calib Thuan Gunn
        if (x2(1).ge.wavbd(1) .and. x2(jmax(i)).le.wavbd(imaxbd)) then
          call sampling(wavbd,flbd,imaxbd,x2,y2,imax2,ls,fs,ts,smax)  
          tmax=0
          do j = 1,smax 
             if (ts(j).ge.tmax) tmax=ts(j)
          enddo   
          do j = 1,smax - 1
            trans     = (ts(j) + ts(j+1))/(2.*tmax)
            flux      = (fs(j)+fs(j+1))/2.
            dlbd      =  ls(j+1)-ls(j) 
            fmeltg(i) = fmeltg(i) + flux*trans*dlbd  
          enddo
c  mtg(*) = mvega(*) + TGcor : TGcor=2.5log(Int[F(BD)TdL]/Int[F(Vega)TdL]) +9.50-0.03
          tg(i) = +2.5*dlog10(fmeltg(i)/fmel(i)) +9.50 -0.03
        else
           tg(i) = -99.99
        endif   
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c Abs mag for sun 
        if (x2(1).ge.wavsun(1).and.x2(jmax(i)).le.wavsun(imaxsun))then
        call sampling(wavsun,lsun,imaxsun,x2,y2,imax2,ls,fs,ts,smax)   
        tmax=0
        arean(i)=0
        fmel(i)=0
        do j = 1,smax 
          if (ts(j).gt.tmax) tmax=ts(j)
        enddo   
        do j = 1,smax - 1
          lmean   = (ls(j) + ls(j+1))/2.
          trans   = (ts(j) + ts(j+1))/(2.*tmax)
          flux    = (fs(j)+fs(j+1))/2.
          dlbd    =  ls(j+1)-ls(j) 
          conv    = c/lmean**2        
          arean(i)= arean(i)+ trans*dlbd*conv
          fmel(i) = fmel(i) + flux*trans*dlbd  
        enddo
c      mo-Mo= 5log D -5 , with D=1U.A. expressed in pc :  mo-Mo= -31.572
c      fo = Lo / 4Pi A^2   with A=1U.A. express in cm
c      mo = -2.5 lg(Lo,nu / 4Pi A^2) -48.59 
c          write(*,*) i,fmel(i),arean(i)
          msun(i)   = -2.5*dlog10(fmel(i)/arean(i)) +68.6227 -48.59
          msun(i)   = msun(i) + 31.572
        else
          msun(i) = -99.99
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  calibrations 
        write(6,30) name(i),id(i),lmoy(i)/10000.,leff(i)/10000.,
     >width(i)/10000.,ab(i),tg(i),veg(i),msun(i),
     > calib(i),l_eff(i)/10000.,fcorr(i)
              
      enddo
 30   format(A10,1x,I3,3(1x,F14.4),4(1x,f7.3),1x,I3,1x,F14.4,1x,f6.3)
c 40   format(A10,2x,f8.3,2x,f8.3,2x,f9.3)
c
      return
 56   write (6,*) 'File ',filtfile2(1:lnblnk(filtfile2)),
     > ' not found -> STOP '
      stop

      end
c
