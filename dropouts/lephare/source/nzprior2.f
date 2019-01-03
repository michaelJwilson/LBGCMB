c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Use prior info on N(z) as a function of magnitude (Iband)
c     and of the type. Work only for color B-I
c
      REAL*8 FUNCTION  nzprior2(Iab,colorRF,z)
      implicit none
c 
      integer*4 mod
      real*8    colorRF
      real*8    z,Iab,pz(3)
      real*8    zmax(3),alpt0(3),ktf(3),ft(3)
      real*8    kt(3),zot(3),pcal,rapp,gammln
      integer*4 choice

      external gammln
c

c     choice=0 for benitez, choice=1 for tolerant fit
      choice=1


      nzprior2 = 1.
      mod     = 1
      if (Iab.le.20) then
         if (z.lt.1) then
              nzprior2 = 1
         else
              nzprior2 = 0
         endif 
         goto 10
      endif
c E/S0
      if (colorRF.ge.1.285) then
      mod     = 1

        if(choice.eq.0)then
          zot(1)  = 0.431
          kt(1)   = 0.091
          alpt0(1)= 2.46
        else
          zot(1)  = 0.431
          kt(1)   = 0.091
          alpt0(1)= 2.46
        endif
      
c Sp
      elseif (colorRF.ge.0.945.and.colorRF.lt.1.285) then
      mod     = 2

        if(choice.eq.0)then
          zot(2)  =0.39 
          kt(2)   = 0.0636
          alpt0(2)= 1.81
        else
          zot(2)  =0.39 
          kt(2)   = 0.1
          alpt0(2)= 1.81
        endif

c Irr
      elseif (colorRF.lt.0.945) then
      mod     = 3

        if(choice.eq.0)then
          zot(3)  = 0.063
          kt(3)   = 0.123
          alpt0(3)=0.91
        else
          zot(3)  = 0.3
          kt(3)   = 0.15
          alpt0(3)=2
        endif

      endif   


c     P(z|T,m0)
      zmax(mod)= zot(mod) + kt(mod)*(Iab-20)
      pz(mod)=z**alpt0(mod)*dexp(-(z/zmax(mod))**alpt0(mod))


c     P(T|m0)
      if(choice.eq.0)then
       ktf(1)=0.147
       ft(1)=0.350
       ktf(2)=0.45
       ft(2)=0.50 
      else
       ktf(1)=0.40
       ft(1)=0.30
       ktf(2)=0.3
       ft(2)=0.35 
      endif

      if(mod.ne.3)then
       rapp=ft(mod)*exp(-ktf(mod)*(Iab-20.d0))
      else
       rapp=1.d0-ft(1)*exp(-ktf(1)*(Iab-20.d0))
     .       -ft(2)*exp(-ktf(2)*(Iab-20.d0))
      endif

      pcal=0.9
c     Normalisation of the probability function
      if(choice.eq.0)then
c       normalization for alpt0(1)= 2.46 alpt0(2)= 1.81 alpt0(3)=0.91
        if(mod.eq.1)pcal=0.8869
        if(mod.eq.2)pcal=0.8891
        if(mod.eq.3)pcal=1.0459
      else
c       normalization for alpt0(1)= 2.46 alpt0(2)= 1.81 alpt0(3)=1.9
        if(mod.eq.1)pcal=0.8869
        if(mod.eq.2)pcal=0.8891
        if(mod.eq.3)pcal=0.8874
      endif
c      pcal=exp(gammln(1.d0/alpt0(mod)+1.d0))
      pcal=zmax(mod)**(alpt0(mod)+1)/alpt0(mod)*pcal

      
c
      nzprior2 = pz(mod) / pcal *rapp
 10   continue

      RETURN
      END
