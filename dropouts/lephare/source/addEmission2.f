      subroutine addEmission2(imagm,Uabs,z,ebv,
     .    DMz,lambf,repf,jmax,ext_em,addem)
c    S.A. : Modified version wrt O.I. version.
c   *  return Emission lines flux [f_nu : in erg/s/cm2/Hz] 
c     in different filters  for a galaxy with absolute magnitude
c     Uabs (not dust corrected) 
c   1: combine all the Emission lines in one spectrum and
c      integrate the spectrum through all filters.
c   2: coefficients between Emission lines have been changed 
c
      implicit none 
c  input
      integer*4  imagm,wmax,maxsize,nbf
c      parameter (wmax=8000)
      parameter  (maxsize=110000)
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_wave.decl'
      integer*4  jmax(nbf)
      real*8     lambf(nbf,maxsize),repf(nbf,maxsize)
      integer*4  nlrep
      real*8     lamb(maxsize),rep(maxsize)
      real*8     Uabs,z,ebv,DMz,c,ext_em(7)
c temp
      integer*4  i,j,is,k,ismax(7),smax
      real*8    ls(maxsize),fs(maxsize),ts(maxsize)
      real*8    fmel,arean,lmean,conv,dlbd,fluxmean,trans
      real*8     lb_line(7),fac_line(7),lmin(7),lmax(7)
      real*8     lbspec(7,wmax),fspec(7,wmax)
      real*8     sigma,lb,lb0,fluxem,fOIInoEx
      parameter (c=2.99792458e+18) 
c   output 
      integer*4  nem
      real*8     lem(wmax),fem(wmax),addem(nbf)
c
ccccccccccccccccccccccccccccccccccc
c
      do k = 1,imagm      
        addEm(k)=0.d0
      enddo
c
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
      fac_line(4)  =  0.61 ! H b  : Moustakas 06: Hb/Ha=0.35->Hb/OII=0.28
      fac_line(5)  =  0.13 ! OIIIa: Moustakas 06: lg[OIII/OII]=-0.33 
      fac_line(6)  =  0.36 ! OIIIb:  ??
      fac_line(7)  =  1.77 ! Ha   : Mouchine 05->STRANGE VALUE!must be 1.77
c Uabs corrected from dust 
      Uabs=Uabs -ebv*ext_em(1)
      fOIInoEx = -0.4*(Uabs+DMz) -6.445
      fOIInoEx = (10.**fOIInoEx)
c
c  distribute the flux through an effective sigma width (10A) 
      sigma=10 
      do j=2,7
         lb0=lb_line(j)*(1+z) 
          fluxem=fOIInoEx*fac_line(j)*10**(-0.4*ebv*ext_em(j))
         lmin(j)=DNINT(lb0)-50.d0
         lmax(j)=DNINT(lb0)+50.d0
         is=0
         do i = 0,100,1
            is=is+1
            lb = DNINT(lb0)-50+i
            lbspec(j,is)=lb
            if (is.eq.1) then 
              fspec(j,is) =0.d0
            else
              fspec(j,is) = 1/(sigma*sqrt(2.*3.14159265))* 
     >                 exp(-(lb-lb0)**2./(2*sigma**2.))*fluxem
c
c      converting in fnu
c              fspec(j,is) = fspec(j,is) *lbspec(j,is)**2 / c
c
              if (fspec(j,is).lt.1.d-60)  fspec(j,is) = 0.d0
            endif
         enddo
         fspec(j,is) =0.d0
         ismax(j) = is 
      enddo
c  Build the spectrum with all Emission lines.
c   Assume that lines are ordered with incr. Lbda
      nem = 1
      lem(1) = 100.d0
      fem(1) = 0.d0
      do j = 2, 7    
         if (j.eq.2) then
            do i = 1,ismax(j)
               nem=nem+1
               lem(nem)=lbspec(j,i)
               fem(nem)=fspec(j,i)
            enddo 
         else
            do i = 1,ismax(j)            
              if (lbspec(j,i).ge.lmin(j-1) .AND.
     >            lbspec(j,i).le.lmax(j-1)) then
                 do k=1,nem
                   if ( ABS(lbspec(j,i)-lem(k)).lt.0.1) then
                      fem(k)=fem(k) + fspec(j,i)
                   endif
                 enddo
              elseif (lbspec(j,i).gt.lmax(j-1)) then
                 nem=nem+1
                 lem(nem)=lbspec(j,i)
                 fem(nem)=fspec(j,i)
              endif
            enddo
        endif 
      enddo 
      nem=nem+1
      lem(nem) = 10000000.d0
      fem(nem) = 0.d0
      if (nem.gt.2000) write(*,*) 'Pb of dimension for EM spectrum'
c
c  compute the AB magnitudes in the filters 
c           
      do i =  1, imagm
         nlrep=jmax(i)
         do j = 1,nlrep
	   lamb(j) = lambf(i,j)
           rep(j) = repf(i,j)
 	 enddo
c    check if lambda start pass fully through the filter
         if (lamb(nlrep).le.lem(nem)) then
            call sampling(lem,fem,nem,lamb,rep,nlrep,ls,fs,ts,smax)
            fmel  = 0.d0
            arean=0.d0
            do j = 1,smax-1
              lmean = (ls(j) + ls(j+1))/2.        ! Average lambda
              conv= c/lmean**2                    ! Conversion in Fnu
              trans = (ts(j) + ts(j+1))/2.        ! Average transmission
              dlbd  =  ls(j+1)-ls(j) 
              arean= arean + trans*dlbd*conv      ! Integrate the filter
              fluxmean = (fs(j)+fs(j+1))/2.       ! average flux
              fmel  = fmel + fluxmean*trans*dlbd  ! Integrate the lfux
            enddo
            addem(i) = fmel/arean                 ! flux in erg/s/cm^2/Hz
          else
            addem(i)   = 0.d0
          endif 
c   add lower limit value  mag(filt)>35 fline[erg/s/cm2/Hz]>10^-33.5
          if (addem(k)<1.e-35) addem(k)=0.             
      enddo
c
      return
      end


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
