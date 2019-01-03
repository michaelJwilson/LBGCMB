c
      subroutine scale_ext(wave,flux,iw,extii,iext,ebv,fext)
c   formatting extinction curve with 
c   respect to the wavelength of SED model ...
      implicit none
c
      integer*4   i,iext
      integer*4   iw,wmax
c      parameter  (wmax=8000)
      INCLUDE 'dim_wave.decl'      
      real*8  wave(wmax),flux(wmax),fout(wmax),ebv,fext(wmax)
      real*8  extil(wmax),extin(wmax),exti(2,wmax)
      real*8  extii(2,wmax)
c
      do i = 1, iext
         extil(i) = extii(1,i)
         extin(i) = extii(2,i)
      end do
      call lbref(wave,iw,extil,extin,iext,fout)
      do i = 1 , iw
          exti(1,i) = wave(i)
          if (wave(i).lt.extil(1) .or.
     >        wave(i).gt.extil(iext)) then
             exti(2,i) = 0.
          else          
             exti(2,i) = fout(i)
          endif
          fext(i) = flux(i)*10**(-0.4*ebv*exti(2,i))
      enddo
      return
      end
c
