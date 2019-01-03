cccccccccccccccccccccccccccccccccccccccccccccc
      subroutine addEmission(imagm,Uabs,z,exti,
     .    DMz,flmoy,flwidth,lambf,repf,jmax,ext_em,addEm)
c
c     Author: Olivier Ilbert (08/02/08)
c     Goal:   Add emission lines according to the absolute magnitude in U 
c             and the extinction
      implicit none
      integer*4 nbf,wmax
      INCLUDE 'dim_filt.decl'    
      INCLUDE 'dim_wave.decl'    
      integer*4  jmax(nbf)
      integer*4  k,imagm
      real*8     c,z,exti,DMz
      real*8     flmoy(nbf),flwidth(nbf)
      real*8     lambf(nbf,wmax),repf(nbf,wmax)
c      real*8     NUV_R
c      integer*4  model
c
      parameter(c=2.99792458e+18)      
      integer*4  nl
      real*8     Uabs,fOIInoEx,lbflux,fline,fluxem
      real*8     lb_line(7),fac_line(7),ext_em(7)
      real*8     addEm(nbf)
      external   fluxem

c     Initialize the emission line flux at 0 in each filter
      do k=1,imagm      
        addEm(k)=0.d0
      enddo
c     lambda Lines 
      lb_line(1) = 2300.   ! NUV 
      lb_line(2) = 1216.   ! Lya
      lb_line(3) = 3727.   ! OII
      lb_line(4) = 4861.   ! Hb
      lb_line(5) = 4959.   ! OIIIa
      lb_line(6) = 5007.   ! OIIIb
      lb_line(7) = 6563.   ! Ha
c     conversion factor between OII and other lines
      fac_line(1)  =  1.   ! UV continuum  ~2300A
      fac_line(2)  =  1.   ! Ly alpha :  no rule 
      fac_line(3)  =  1.   ! OII  : Kennicutt 98: UVcont -> OII
      fac_line(4)  =  0.28 ! H b  : Moustakas 06: Hb/Ha=0.35->Hb/OII=0.28
      fac_line(5)  =  0.30 ! OIIIa: Moustakas 06: lg[OIII/OII]=-0.33 
      fac_line(6)  =  0.47 ! OIIIb:  ??
      fac_line(7)  =  0.80 ! Ha   : Mouchine 05->STRANGE VALUE!must be 1.77
c
c     NUV Absolute magnitude dust corrected  -> converted in OII flux
      Uabs     = Uabs-exti*ext_em(1)
      fOIInoEx = -0.4*(Uabs+DMz) -6.345
      fOIInoEx = (10.**fOIInoEx)
c    
      do k=1,imagm              ! loop on filters 
c        if(model.lt.8)  addEm(k) = 0
c        if(model.lt.8)  goto 78  !  Apply it only on blue galaxies
c      
        do nl= 2,7 
          if(dabs(flmoy(k)-(lb_line(nl)*(1.+z))).le.flwidth(k))then
             lbflux=lb_line(nl)*(1.+z)
             fline=fOIInoEx*fac_line(nl)*10**(-0.4*exti*ext_em(nl))
c            Integrate the line flux in the filter 
             addEm(k)=fluxem(lbflux,fline,k,lambf,repf,jmax)
c             write(*,*) k,nl,z,exti,fOIInoEx,fline,addEm(k)
c   add lower limit value  mag(filt)>35 fline[erg/s/cm2/Hz]>10^-35
             if (addEm(k)<1.d-50) addEm(k)=0.             
          endif
        enddo
c 78     continue
      enddo
c
      return
      end
c
cc NOTE from SA  based on Kennicutt 98, ARAA, 36   for a Salpeter IMF, dust corrected
c   SFR=1.4 10^{-28}  Lnu  [erg/s/Hz]  = 7.9 10^{-42} L(Ha) [erg/s] = 1.4 10^{-41} L(OII) [erg/s]   
c    
c  By integrating L[OII] over 10A  : L(OII)[erg/s] = L(OII)[erg/s/Hz] *Dnu = L(OII)*c/lbda^2 * Dlbda 
c  assuming :lbda=3727, Dlbda=10   : L(OII)[erg/s] = 2.16 10^12 L(OII)[erg/s/Hz] 
c  Replacing in Kennicutt formulae : 
c                    1.4 10^{-28}  Lnu  [erg/s/Hz] =  1.4 10^{-41}  2.16 10^12 L(OII)[erg/s/Hz]
c
c   L(OII) [erg/s/Hz] = 4.63 Lnu(UV) [erg/s/Hz]   -> log(f(OII)/f(nuv)) = 0.66
c
c  lg(fnuv)[erg/s/cm2/Hz] = -0.4 (NUVabs + DM(z) + 48.59) = lg(f(OII)) -0.66  
c  lg(fOII)[erg/s/cm2/Hz] = -0.4(NUVabs + DM(z)) -18.78 
c  lg(fOII))[erg/s/cm2]   = lg(2.16 10^12) + lg(fOII)[erg/s/cm2/Hz] 
c  lg(fOII))[erg/s/cm2]   = lg(2.16 10^12) -0.4(NUVabs + DM(z)) -18.78 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  The conversion absolute magnitude in NUV to OII flux line is : c
c  lg(fOII))[erg/s/cm2]   = -0.4(NUVabs + DM(z)) -6.445           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              OTHER Emission Lines  :
c  * For Halpha with Eq above, it is strightforward :
c   -->   f(Ha)/f(OII)= 1.77   !   you assume 0.8 
c
c  * For Hbeta  based on the theoretical ratio: f(Ha)/f(Hb)=2.91 (McCall 1982)
c   -->   f(Hb)/f(OII) = 0.61  !  you assume  0.28
c   
c  * For OIII : In the plot AGN/SF   based on lg(OIII[5007]/Hb) vs lg(NII/Ha) 
c    for SF galaxies  OIII/Hb varies by a factor 30 !!
c    while NII/Ha coulb be considered around NII/Ha=0.32 +/-0.2
c    Assuming OIII/Hb in the middle ->  lg(OIII/Hb)=-0.2 -> OIII[5007]/Hb=0.6
c
c   -->   f(OIII[5007])/f(OII)=0.36  ! you assume 0.47 
c
c    With McCall (85) relation : OIII[4959+5007]/Hb =1.35 OIII[5007]/Hb
c
c   -->   f(OIII[4959])/f(OII) = 0.13  ! you assume 0.30
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  with extinction Ha/Hb)ext = Ha/Hb)th 10^(-0.4*EB-V*(Aa-Ab))
c
