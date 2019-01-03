c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Use prior info on N(z) as a function of magnitude (Iband)
c
      REAL*8 FUNCTION  nzprior3(Iab,z)
      implicit none
c
      real*8    z,Iab,pz
      real*8    zmax,alpt0
      real*8    kt,zot,pcal,gammln

      external gammln
c

c     No prior if Iab<20
      nzprior3 = 1.
      if (Iab.le.20) then
         if (z.lt.1) then
              nzprior3 = 1
         else
              nzprior3 = 0
         endif 
         goto 10
      endif

c     parameterization of the z distribution (functional form of Benitez)
      zot  = 0.3
c original value from Ilbert 
c      kt   = 0.14
c      alpt0= 1.2
c  slightly modified by myself 
      kt   = 0.1
      alpt0= 1.4
c     P(z|m0)
      zmax= zot + kt*(Iab-20)
      pz=z**alpt0*dexp(-(z/zmax)**alpt0)

c     Normalization (always same gamma)
c      pcal=exp(gammln(1.d0/alpt0+1.d0))
c      pcal = 0.887363315891821
      pcal=0.911423338  
      pcal=zmax**(alpt0+1)/alpt0*pcal

      
c     prior to apply
      nzprior3 = pz / pcal 

 10   continue

      RETURN
      END
