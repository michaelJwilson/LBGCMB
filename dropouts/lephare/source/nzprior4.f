c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Use prior info on N(z) as a function of magnitude (Iband)
c     and of the type. Work only for color B-I
c
      REAL*8 FUNCTION  nzprior4(Iab,modele,z)
      implicit none
c
      integer*4 mod,modele
      real*8    z,Iab,pz(4)
      real*8    zmax(4),alpt0(4),ktf(4),ft(4)
      real*8    kt(4),zot(4),pcal,rapp,gammln

      external gammln
c

      nzprior4 = 1.
      mod     = 1
      if (Iab.le.20) then
         if (z.lt.1) then
              nzprior4 = 1
         else
              nzprior4 = 0
         endif 
         goto 10
      endif
c E/S0
      if (modele.le.15) then
         mod     = 1
         zot(1)  = 0.42113 
         kt(1)   = 0.14364
         alpt0(1)= 2.19517

c Sbc
      elseif (modele.gt.15..and.modele.le.30) then
         mod     = 2
         zot(2)  = 0.34897
         kt(2)   = 0.12844
         alpt0(2)= 2.11268 



c Scd
      elseif (modele.gt.30.and.modele.le.45) then
         mod     = 3
         zot(3)  = 0.24302
         kt(3)   = 0.12937 
         alpt0(3)= 1.64212

c Irr
      elseif (modele.gt.45.and.modele.le.100) then
         mod     = 4
         zot(4)  = 0.18636 
         kt(4)   = 0.14101
         alpt0(4)= 1.29755



      endif   

c     P(z|T,m0)
      zmax(mod)= zot(mod) + kt(mod)*(Iab-20)
      pz(mod)=z**alpt0(mod)*dexp(-(z/zmax(mod))**alpt0(mod))


c     P(T|m0)
      ktf(1)=   0.50232
      ft(1) =   0.45110
      ktf(2)=   0.41061
      ft(2) =   0.19304
      ktf(3)=   0.15240
      ft(3) =   0.33146 
      ktf(4)=  -0.43913 
      ft(4) =   0.15093


      if(mod.ne.4)then
       rapp=ft(mod)*exp(-ktf(mod)*(Iab-20.d0))
      else
       rapp=1.d0-ft(1)*exp(-ktf(1)*(Iab-20.d0))
     .          -ft(2)*exp(-ktf(2)*(Iab-20.d0))
     .          -ft(3)*exp(-ktf(3)*(Iab-20.d0))
      endif


c     Normalisation of the probability function
c      pcal=exp(gammln(1.d0/alpt0(mod)+1.d0))
      pcal = 0.89 
      if(mod.eq.1)pcal=   0.88562 
      if(mod.eq.2)pcal=   0.88566 
      if(mod.eq.3)pcal=   0.89456
      if(mod.eq.4)pcal=   0.92393
      pcal=zmax(mod)**(alpt0(mod)+1)/alpt0(mod)*pcal


c
      nzprior4 = pz(mod) / pcal *rapp
 10   continue

      RETURN
      END
