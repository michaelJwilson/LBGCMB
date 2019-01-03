c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Use prior info on N(z) as a function of magnitude (Iband)
c     and of the type. Work only for color B-I
c
      REAL*8 FUNCTION  nzpriorVVDS5(Iab,modele,z)
      implicit none
c
      integer*4 mod,modele
      real*8    z,Iab,pz(4)
      real*8    zmax(4),alpt0(4),ktf(4),ft(4)
      real*8    kt(4),zot(4),pcal,rapp,gammln

      external gammln
c

      nzpriorVVDS5 = 1.
      mod     = 1
      if (Iab.le.20) then
         if (z.lt.1) then
              nzpriorVVDS5 = 1
         else
              nzpriorVVDS5 = 0
         endif 
         goto 10
      endif
c E/S0
      if (modele.le.15) then
         mod     = 1
         zot(1)  = 0.45181 
         kt(1)   = 0.13677
         alpt0(1)= 3.33078 

c Sbc
      elseif (modele.gt.15..and.modele.le.30) then
         mod     = 2
         zot(2)  = 0.16560
         kt(2)   = 0.12983 
         alpt0(2)= 1.42815



c Scd
      elseif (modele.gt.30.and.modele.le.45) then
         mod     = 3
         zot(3)  = 0.21072
         kt(3)   = 0.14008
         alpt0(3)= 1.58310

c Irr
      elseif (modele.gt.45.and.modele.le.100) then
         mod     = 4
         zot(4)  = 0.20418
         kt(4)   = 0.13773 
         alpt0(4)= 1.34500

      endif   

c     P(z|T,m0)
      zmax(mod)= zot(mod) + kt(mod)*(Iab-20)
      pz(mod)=z**alpt0(mod)*dexp(-(z/zmax(mod))**alpt0(mod))


c     P(T|m0)
      ktf(1)=    0.47165  
      ft(1) =    0.43199
      ktf(2)=    0.30663
      ft(2) =    0.07995 
      ktf(3)=    0.12715
      ft(3) =    0.31162
      ktf(4)=    -0.34437
      ft(4) =    0.21220


      if(mod.ne.4)then
       rapp=ft(mod)*exp(-ktf(mod)*(Iab-20.d0))
      else
       rapp=1.d0-ft(1)*exp(-ktf(1)*(Iab-20.d0))
     .          -ft(2)*exp(-ktf(2)*(Iab-20.d0))
     .          -ft(3)*exp(-ktf(3)*(Iab-20.d0))
       if(rapp.le.0)then
        write(*,*)"PROBLEME PRIOR, SUM>1",rapp
        rapp=0.d0
       endif
      endif


c     Normalisation of the probability function
c      pcal=exp(gammln(1.d0/alpt0(mod)+1.d0))
      if(mod.eq.1)pcal=  0.89744  
      if(mod.eq.2)pcal=  0.90868  
      if(mod.eq.3)pcal=  0.89747
      if(mod.eq.4)pcal=  0.91760 
      pcal=zmax(mod)**(alpt0(mod)+1)/alpt0(mod)*pcal


c
      nzpriorVVDS5 = pz(mod) / pcal *rapp
 10   continue

      RETURN
      END
