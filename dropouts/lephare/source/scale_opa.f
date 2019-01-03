c
      subroutine scale_opa(wav,fl,iw,opal,opat,nopa,zobs,
     >                     wnew,fnew,kmax)
c  Apply the opacity at a given zobs 
c  to the SED fluxes by resampling both curves
c  according to their wavelengthes 

      implicit none
      integer*4 i,n,k,kmax,iw,nop,nopa(81),nmax,wmax
c      parameter (wmax=8000)
      INCLUDE 'dim_wave.decl'      
      real*8    val,zobs
      real*8    wav(wmax),fl(wmax),wnew(wmax),fnew(wmax)
      real*8    opa(wmax),wopa(wmax),opal(81,wmax),opat(81,wmax)
      real*8    wn(wmax),fn(wmax),on(wmax)
c     
            val = zobs/0.1 + 1
            nop = idnint(val)
c            nop= 61
c            write(*,*) zobs,val,nop
            if (nop.gt.81) nop = 81
            k = 0
            do n = 1, nopa(nop), 3    ! I use double lbda step (Dlbda=2A)
c              if ((zobs/0.1+1) .lt. val) then
c                 if (opat(nop,n).le.1.e-5)   opat(nop,n)=0.
c                 if (opat(nop+1,n).le.1.e-5) opat(nop+1,n)=0.
c                  k = k + 1
c                  opa(k)= (opat(nop,n) + opat(nop+1,n))/2.
c                  wopa(k) = opal(nop,n)
c              else 
                  if (opat(nop,n).le.1.e-5)   opat(nop,n)=0.
                  k = k + 1 
                  opa(k)= opat(nop,n)
                  wopa(k) = opal(nop,n)
c              endif
           enddo  
           kmax = k  
c  resampling the 20->1215A range
           call sampling(wav,fl,iw,wopa,opa,kmax,wn,fn,on,nmax)
c  re-defined the flux of fl and wave with the new wn and fn*on 
           k=0
           do i = 1,nmax
              k = k + 1
              wnew(k) = wn(i)
c              fnew(k) = fn(i)
              fnew(k) = fn(i)*on(i)
           enddo
           do i = 1, iw
              if (wav(i).gt.wn(nmax)) then
                 k = k + 1
                 wnew(k) = wav(i)
                 fnew(k) = fl(i)
              endif
           enddo   
           kmax=k
c
      return
      end
c
