c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Use prior info on N(z) as a function of magnitude (Iband)
c     and of the type. WORK only for CWW libraries (4-6 spectra)
c
      REAL*8 FUNCTION  nzprior(Iab,model,z)
      implicit none
c
      integer*4 mod,model
      real*8    z,Iab,pz(3),ptyp(3)
      real*8    zmax(3),alpt0(3),alpt1(3)
      real*8    kt(3),zot(3),zc,pcal
c
      nzprior = 1.
      mod     = 1
      if (Iab.le.20) then
         if (z.lt.1) then
              nzprior = 1
         else
              nzprior = 0
         endif 
         goto 10
      endif
c E/S0
      if (model.le.2) then
      mod     = 1
      zot(1)  = 0.48
      kt(1)   = 0.061
      alpt0(1)= 2.26 
c      alpt1(1) = 2.26
      alpt1(1)= 1.8
      ptyp(1) = 0.35*dexp(-0.47*(Iab-20))
c Sp
      elseif (model.gt.2.and.model.le.4) then
      mod     = 2
      zot(2)  = 0.44
      kt(2)   = 0.044
      alpt0(2)= 1.71
c      alpt1(2)= 1.71
      alpt1(2)= 1.11
      ptyp(2) = 0.5*exp(-0.165*(Iab-20)) 
c Irr
      elseif (model.gt.4) then
      mod     = 3
      zot(3)  = 0.038
      kt(3)   = 0.178
      alpt0(3)= 1.125
c      alpt1(3)= 0.75
      alpt1(3)= 0.7
      ptyp(3) = 1 - ptyp(1) -ptyp(2)
      endif   
      zmax(mod)= zot(mod) + kt(mod)*(Iab-20)
      pz(mod)=z**alpt0(mod)*dexp(-(z/zmax(mod))**alpt1(mod))
c
      zc = (alpt0(mod)/alpt1(mod))**(1/alpt1(mod)) * zmax(mod)
      pcal = zc**alpt0(mod) * dexp(-(zc/zmax(mod))**alpt1(mod))
      nzprior = pz(mod) / pcal
 10   continue
      RETURN
      END
