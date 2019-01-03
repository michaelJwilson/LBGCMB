      real*8 function fluxem(lbflux,fluxin,filtre,
     .          lambf,repf,jmax)

      implicit none      
      integer*4  nbf,maxsize,wmax
      real*8     c
      parameter(c=2.99792458e+18)      
      parameter  (maxsize=110000)
      include 'dim_filt.decl'
      include 'dim_wave.decl'
c      real*8    lambf(nbf,maxsize),repf(nbf,maxsize),sigma,sum
c      real*8    lamb(maxsize),rep(maxsize)
c      real*8    lbem(maxsize),flem(maxsize)
c      real*8    ls(maxsize),fs(maxsize),ts(maxsize)
      real*8    lambf(nbf,wmax),repf(nbf,wmax),sigma,sum
      real*8    lamb(wmax),rep(wmax)
      real*8    lbem(wmax),flem(wmax)
      real*8    ls(wmax),fs(wmax),ts(wmax)
      integer*4 nlrep,j,smax
      integer*4 jmax(nbf),filtre,nlf    
      real*8    fmel,arean,lmean,trans,fluxin,fluxmean
      real*8    conv,dlbd,lbflux
c
c     number of lines, lambda, transmission of the considered filter
      nlrep=jmax(filtre) 
      do j = 1,nlrep
          lamb(j) = lambf(filtre,j)
          rep(j)  = repf(filtre,j)
      enddo
c     create a dirac emission line on a grid of lambda every 1 A 
c      and limited by the filter lambda range
      sigma=10
      sum=0.d0
      nlf=lamb(nlrep)-lamb(1)
      do j=1,nlf
         lbem(j)=lamb(1)+j-1
c        distribute the flux according to a gaussian with sigma A
         flem(j)=1/(sigma*sqrt(2.*3.14159265)) * 
     .            exp(-(lbem(j)-lbflux)**2./(2*sigma**2.))*fluxin
         sum=sum+flem(j)
      enddo
c     resample with the 1A grid
      call sampling(lbem,flem,nlf,lamb,rep,nlrep,ls,fs,ts,smax)
      fmel  = 0.d0
      arean=0.d0
      do j = 1,smax-1
          lmean = (ls(j) + ls(j+1))/2.         ! Average lambda
          conv= c/lmean**2                     ! Conversion in Fnu
          trans = (ts(j) + ts(j+1))/2.         ! Average transmission
          dlbd  =  ls(j+1)-ls(j) 
          arean= arean + trans*dlbd*conv       ! Integrate the filter
          fluxmean = (fs(j)+fs(j+1))/2.        ! average flux
          fmel  = fmel + fluxmean*trans*dlbd   ! Integrate the lfux
      enddo
      fluxem  = fmel/arean                     ! flux in line in erg/s/cm^2
c
      return
      end      
