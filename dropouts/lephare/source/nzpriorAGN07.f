c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Use prior info on N(z)
c     Only valid for AGN Type 1
c

      REAL*8 FUNCTION  nzpriorAGN07(z,agn_a,agn_zm)
      implicit none
c
      real*8    z,agn_a,agn_zm
      real*8    pz,pcal
      real*8    gammln

      external gammln
c

      nzpriorAGN07 = 1.


      pz=z**agn_a*dexp(-(z/agn_zm)**agn_a)

c     Normalisation of the probability function
c      pcal=exp(gammln(1.d0/agn_a+1.d0))
c      write(*,*)"PCAL",pcal
      pcal= 0.892244502
      pcal=agn_zm**(agn_a+1)/agn_a*pcal

c
      nzpriorAGN07 = pz / pcal 

c 10   continue

      RETURN
      END
