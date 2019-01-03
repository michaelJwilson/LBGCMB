c   last modif : 06/07/00
      SUBROUTINE read_spec(lib,rec,extil,ebv,zb,
     >exti,iext,opal,opat,iopa,
     >wavemin,wavemax,sampl,wav,fl,spmax)
c
c     this program extracts the spectrum for QSO and Galaxy 
c     from the library and put it at the obs. z and used extinction ...
c     It return flux in frequency fnu(z) and wave redshifted wave(z)
      implicit none
      integer*4  k,j,inmax,rec
      integer*4  test,wmax,spmax
      integer*4  reclmax,extil
      integer*4  lnblnk,nr,nmod,iw,sampl
      integer*4  iopa(81),iext(10),iextinc,libfmt
c      
c      parameter  (wmax=8000)
      INCLUDE 'dim_wave.decl'      
c
      real*8     exti(10,2,wmax),opal(81,wmax),opat(81,wmax)
      real*8     extinc(2,wmax),extis(2,wmax)
      real*8     wave(wmax),wn(wmax),fn(wmax)
      real*8     dummy,vec(wmax),wav(wmax),fl(wmax)
      real*8     c,ebv,zb,wavemax,wavemin
      parameter (c=2.99792458e+18) 
      character*4096  file,lib,zpdir,zpwork,param,paravc(500)
c   ... INPUT ...
      call getenv('LEPHAREDIR',zpdir)
      test=lnblnk(zpdir)
      if (test .eq. 0) then
        write(*,*) 'WARNING :  variable LEPHAREDIR not defined'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(6,*) 'WARNING :  variable LEPHAREWORK not defined'
        stop
      endif
c
      do k = 1, wmax 
         wav(k)= 0
         fl(k) = 0
      enddo
      file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     > // lib(1:lnblnk(lib)) // '.doc'
      param='RECORD_LENGTH'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') reclmax
c
      libfmt=0
      param='LIB_FMT'
      call read_para(param,file,paravc,test)
      if (test.eq.1 .AND.
     >     paravc(1)(1:lnblnk(paravc(1))).eq.'BC_STOCH') libfmt=1
      
c     
      file = zpwork(1:lnblnk(zpwork)) // '/lib_bin/'
     >       // lib(1:lnblnk(lib)) // '.bin'
      open(28,file=file,status='unknown',access='direct',recl=reclmax)
c pick the spectrum 
      read(28,rec=rec) nr,nmod,dummy,iw,(vec(k),k=1,iw)
c      write(6,*) nr,nmod,dummy,iw,lib(1:lnblnk(lib))
c search the wavelengths :
      if (libfmt .eq. 0 ) then 
        do j = rec-1, 1,-1
           read(28,rec=j) nr,nmod,dummy,iw,(wave(k),k=1,iw)
           if (dummy.eq.-1.0)   goto 33
        enddo
      else
c wavelengths in first row for BC_STOCH Library 
        read(28,rec=1) nr,nmod,dummy,iw,(wave(k),k=1,iw)
      endif
 33   close(28)
c   add  IGM opacity
      call scale_opa(wave,vec,iw,opal,opat,iopa,zb,wn,fn,inmax)
c add extinction 
      if (ebv.gt.0) then
         iextinc=iext(extil)
         do k = 1, iextinc
            extinc(1,k)=exti(extil,1,k)
            extinc(2,k)=exti(extil,2,k)
         enddo   
         call lambda(wn,inmax,extinc,iextinc,extis)
       endif  
       spmax=0
       do k = 1, inmax ,sampl
          if ( (wn(k)*(1+zb)).le.wavemax .and.
     >         (wn(k)*(1+zb)).ge.wavemin )  then 
             spmax = spmax + 1
             wav(spmax)=wn(k)*(1+zb)
             fn(k)=fn(k)/(1+zb)
             if (ebv.gt.0)  fn(k)=fn(k)*10**(-0.4*ebv*extis(2,k))
             if (fn(k).le.0) then
                fl(spmax)= 150.
             else
             fl(spmax)=-2.5*dlog10(fn(k)*wav(spmax)*wav(spmax)/c)-48.59
             endif
          endif
      enddo
c
      RETURN
      end
c
