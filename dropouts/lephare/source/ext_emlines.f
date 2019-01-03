c
cccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ext_emlines(exti,iext,lUV,fUV,jUV,ext_em)
c
c     Author: Olivier Ilbert (08/02/08)
c     Goal:   compute extinction according to the different extinction curve
c             in NUV and in the different emission lines

      integer*4  wmax
c      parameter (wmax=8000)
      INCLUDE 'dim_wave.decl'      
c
      integer*4 jUV,iext,i,k
c      integer*4 ifilt
      real*8 lext(wmax),ext(wmax)
      real*8 ext_em(7),wav_em(7),lUV(wmax),fUV(wmax),exti(2,wmax)
      real*8 lfilti(wmax),filti(wmax),aint
c
c     central wavelength of each line
      wav_em(1)=2300      ! NUV 
      wav_em(2)=1216      ! Lya
      wav_em(3)=3727      ! OII 
      wav_em(4)=4861      ! Hb
      wav_em(5)=4959      ! OIII
      wav_em(6)=5007      ! OIII 
      wav_em(7)=6563      ! Ha
c      
c     Extinction in a vector 1 dimension
      do i=1,iext
         lext(i)=exti(1,i)
         ext(i)=exti(2,i)
      enddo
c    
c     Loop on the UV filter and each emission line      
      do k=1,7
c
c       NUV filter
        if(k.eq.1)then
          ifilti=jUV
          do i=1,ifilti
            lfilti(i)=lUV(i)
            filti(i)=fUV(i)
          enddo
        endif
c
c       For each emission line, create a filter
        if(k.gt.1)then
          ifilti=200
          do i=1,ifilti
            lfilti(i)=wav_em(k)+dble(i)-100.
            filti(i)=0.d0
            if(abs(lfilti(i)-wav_em(k)).lt.20) filti(i)=1.
          enddo
        endif
c
c       Compute the extinction
        call comp_ext2(lfilti,filti,ifilti,lext,ext,iext,aint)
c
        ext_em(k)=aint
c        write(*,*)"Extinction ",aint
      enddo
c   
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     From Stephane, filter_extinc.f, slighly modified
c
      subroutine comp_ext2(lfilt,filt,ifilt,lext,ext,iext,aint)
      implicit none 
      integer*4 wmax,ifilt,i,iext,smax
      INCLUDE 'dim_wave.decl'      
c      parameter (wmax=8000)
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
         fint =fint+( ts(i) + ts(i+1) ) /2.  * (ls(i+1) - ls(i))
         aint =aint+(ts(i)+ts(i+1))*(ls(i+1)-ls(i))*(fs(i)+fs(i+1))/4.
      enddo
      aint = aint / fint 
c
      return
      end
