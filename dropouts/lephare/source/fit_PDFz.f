ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     19/09/04
c     Fit the PDFz with 4 Gaussian functions
c     The minimization is done using the minuit library.
      subroutine fit_PDFz(xp,yp,chimax,nb_coef,zb,nb,PDFz_coef)

c     xp yp : data to fit    
c     chimax : number of lines to read in yp
c     nb_coef : number of coefficient for the fit
c     zb nb : z value of the secondary peaks
c     PDFz_coef* : output coefficients
c

      implicit none

      integer*4 chisize,i,k
      INCLUDE 'dim_zchi.decl'

      integer*4 chimax,chiPDFz,nb_coef,ncoef,nb(500),first_pass
      real*8 xp(chisize),yp(chisize),xPDFz(chisize),yPDFz(chisize)
      real*8 PDFz_coef(500),zb(chisize)      

c     variable minuit
      real*8 vstrt(500),stp(500),zero,u(500)
      integer*4 nprm(500),ierflg
      character*4096 pnam(500)
      data zero /0./      

c     Function for the minimization
      external fcnPDFz2

c     common  minuit
      common/func_tabPDZ/xPDFz,yPDFz
      common/func_intPDZ/chiPDFz,ncoef
      common/mn7ext/u

c     Common variable for minuit
      chiPDFz=chimax
      ncoef=nb_coef
      do k=1,chimax
        xPDFz(k)=xp(k) 
        yPDFz(k)=yp(k) 
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Initialize the free parameters : 
c     Try the be the close as possible of the true values
c
c     number, name ,step
      do i=1,nb_coef
        nprm(i)=i
        pnam(i)='a0'//char(i+47)
        stp(i)=0.02d0
      enddo
c     Take the value of the PDFz at each peak
      vstrt(1)=yp(nb(1))
      vstrt(4)=yp(nb(2))
      vstrt(7)=yp(nb(3))
      vstrt(10)=yp(nb(4))
c     If no secondary peak, set the normalization at 0
      if(zb(nb(1)).le.0)vstrt(1)=0.d0
      if(zb(nb(2)).le.0)vstrt(4)=0.d0
      if(zb(nb(3)).le.0)vstrt(7)=0.d0
      if(zb(nb(4)).le.0)vstrt(10)=0.d0
c     Set z parameter to the value of the z secondary peaks
      vstrt(2)=zb(nb(1))
      vstrt(5)=zb(nb(2))
      vstrt(8)=zb(nb(3))
      vstrt(11)=zb(nb(4))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Two iterations : the first one with only the sigma parameters free
c     the second one with all free
      first_pass=0
 10   continue

c     Initialize the value of sigma
      if(first_pass.eq.0)then
c      arbitrary
       vstrt(3)=0.1d0
       vstrt(6)=0.1d0
       vstrt(9)=0.1d0
       vstrt(12)=0.1d0
      else
c      Initialize with the best value after the 1st iteration
       vstrt(3)=U(3)
       vstrt(6)=U(6)
       vstrt(9)=U(9)
       vstrt(12)=U(12)
      endif

c     Initialize the minuit output, write the less as possible
      call mninit(5,6,7)
      call mncomd(fcnPDFz2,'set printout -1',ierflg,0)
      call mncomd(fcnPDFz2,'set nowarnings',ierflg,0)
c     Minuit parameters initialization
      do i= 1,nb_coef
       call mnparm(nprm(i),pnam(i),vstrt(i),stp(i),zero,zero,ierflg)
       if (ierflg .ne. 0)  then
        write (6,'(a,i3)')  ' unable to define parameter no.',i
        stop
       endif
      enddo
c     Best minuit accuracy
      call mncomd(fcnPDFz2,'set str 2',ierflg,0)
c     Limits allowed
      call mncomd(fcnPDFz2,'set lim 1 -1.d-2 10000.',ierflg,0)
      call mncomd(fcnPDFz2,'set lim 3 0.d0 1.',ierflg,0)
      call mncomd(fcnPDFz2,'set lim 4 -1.d-2 10000.',ierflg,0)
      call mncomd(fcnPDFz2,'set lim 6 0.d0 1.',ierflg,0)
      call mncomd(fcnPDFz2,'set lim 7 -1.d-2 10000.',ierflg,0)
      call mncomd(fcnPDFz2,'set lim 9 0.d0 10.',ierflg,0)
      call mncomd(fcnPDFz2,'set lim 10 -1.d-2 10000.',ierflg,0)
      call mncomd(fcnPDFz2,'set lim 12 0.d0 11.',ierflg,0)

cc     To save time : never fit the z 
cc     but problem with the asymmetric distribution
c      call mncomd(fcnPDFz2,'fix 2',ierflg,0)
c      call mncomd(fcnPDFz2,'fix 5',ierflg,0)
c      call mncomd(fcnPDFz2,'fix 8',ierflg,0)


c     Fix the z and normalization values in the first pass
c     leave only free the sigma
      if(first_pass.eq.0)then
       call mncomd(fcnPDFz2,'fix 1',ierflg,0)
       call mncomd(fcnPDFz2,'fix 2',ierflg,0)
       call mncomd(fcnPDFz2,'fix 4',ierflg,0)
       call mncomd(fcnPDFz2,'fix 5',ierflg,0)
       call mncomd(fcnPDFz2,'fix 7',ierflg,0)
       call mncomd(fcnPDFz2,'fix 8',ierflg,0)
      endif

c     If no secondary peak, no fit at all
      if(zb(nb(1)).le.0)then
          call mncomd(fcnPDFz2,'fix 1',ierflg,0)
          call mncomd(fcnPDFz2,'fix 2',ierflg,0)
          call mncomd(fcnPDFz2,'fix 3',ierflg,0)
      endif
      if(zb(nb(2)).le.0)then
          call mncomd(fcnPDFz2,'fix 4',ierflg,0)
          call mncomd(fcnPDFz2,'fix 5',ierflg,0)
          call mncomd(fcnPDFz2,'fix 6',ierflg,0)
      endif
      if(zb(nb(3)).le.0)then
          call mncomd(fcnPDFz2,'fix 7',ierflg,0)
          call mncomd(fcnPDFz2,'fix 8',ierflg,0)
          call mncomd(fcnPDFz2,'fix 9',ierflg,0)
      endif
      if(zb(nb(4)).le.0)then
          call mncomd(fcnPDFz2,'fix 10',ierflg,0)
          call mncomd(fcnPDFz2,'fix 11',ierflg,0)
          call mncomd(fcnPDFz2,'fix 12',ierflg,0)
      endif
c     Minimization
      call mnexcm(fcnPDFz2,'migrad',zero,0,ierflg,0)


c     If the sigma is determined, return at the beginning
c     and fit all free parameters 
      if(first_pass.eq.0)then
         first_pass=1
         goto 10
      endif



c     output values of the fit
      do k=1,nb_coef
       PDFz_coef(k)=u(k)
      enddo

      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     function call by minuit PDFz to compute the chi2 
      subroutine fcnPDFz2(npar,gin,chi,x,iflag)

      implicit none

      integer*4 chisize,i
      INCLUDE 'dim_zchi.decl'

      integer*4 chiPDFz,npar,ncoef
      real*8 xPDFz(chisize),yPDFz(chisize)
      real*8 chi,x(npar),yf

      real*8 gin
      integer*4 iflag

      common/func_tabPDZ/xPDFz,yPDFz
      common/func_intPDZ/chiPDFz,ncoef


c     initialize chi2
      chi=0.d0
c     loop over PDFz
      do i=1,chiPDFz
c       Functional form : sum of 4 Gaussian functions
        yf=x(1)*exp(-((xPDFz(i)-x(2))/x(3))**2.d0/2.)
     .    +x(4)*exp(-((xPDFz(i)-x(5))/x(6))**2.d0/2.)
     .    +x(7)*exp(-((xPDFz(i)-x(8))/x(9))**2.d0/2.)
     .    +x(10)*exp(-((xPDFz(i)-x(11))/x(12))**2.d0/2.)
c        chi2
         chi=chi+(yPDFz(i)-yf)**2.d0
      enddo

      end
