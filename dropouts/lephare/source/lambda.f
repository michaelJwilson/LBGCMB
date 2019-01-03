c
      subroutine lambda(wave,iw,exti,iext,extis)
c   formatting extinction curve with 
c   respect to the wavelength of SED model ...
      implicit none
c
      integer*4  i,iext
      integer*4  iw,wmax
c      parameter  (wmax=8000)
      INCLUDE 'dim_wave.decl'      
      real*8     wave(wmax),fout(wmax)
      real*8     extil(wmax),extin(wmax),extis(2,wmax),exti(2,wmax)
c
      do i = 1, iext
         extil(i) = exti(1,i)
         extin(i) = exti(2,i)
       end do
       call lbref(wave,iw,extil,extin,iext,fout)
       do i = 1 , iw
          extis(1,i) = wave(i)
          if (wave(i).lt.extil(1) .or.
     >        wave(i).gt.extil(iext)) then
             extis(2,i) = 0.
          else          
             extis(2,i) = fout(i)
          endif
       enddo
      return
      end
c
