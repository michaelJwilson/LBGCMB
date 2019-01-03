      subroutine addemlines(Uabs,z,DMz,aext_lb,elaw,ebv,frac,
     >       emspec,nem)
c     Return Emission lines spectrum in fnu [erg/s/cm2/Hz]
      implicit none 
c  input
      integer*4  elaw,wmax
c      parameter (wmax=8000)
      INCLUDE 'dim_wave.decl'      
      real*8     Uabs,z,ebv,aext_lb(10,7),DMz,frac(6)
c temp
      integer*4  i,j,is,k,ismax(7)
      real*8     lb_line(7),fac_line(7),lmin(7),lmax(7)
      real*8     lbspec(7,wmax),fspec(7,wmax)
      real*8     sigma,lb,lb0,fem,fOIInoEx,c
      parameter (c=2.99792458e+18) 
c   output 
      integer*4  nem
      real*8     emspec(2,wmax)
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
c Uabs is dust free lumisosity from model
      fOIInoEx = -0.4*(Uabs+DMz) -6.445
      fOIInoEx = (10.**fOIInoEx)
c
      sigma=10 
      nem=0
      do j=2,7
         lb0=lb_line(j)*(1+z) 
         fem=fOIInoEx*fac_line(j)*10**(-0.4*ebv*aext_lb(elaw,j))
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
     >                 exp(-(lb-lb0)**2./(2*sigma**2.))*fem
              fspec(j,is) = fspec(j,is) *lbspec(j,is)**2 / c
              if (fspec(j,is).lt.1.d-60)  fspec(j,is)= 0.d0
            endif
         enddo
         fspec(j,is) =0.d0
         ismax(j) = is 
      enddo
c  Sum the spectra at same wavelengths. Assume that lines are ordered with incr. Lbda
      nem = 1
      emspec(1,1) = 100.d0
      emspec(2,1) = 0.d0
      do j = 2, 7    
         if (j.eq.2) then
            do i = 1,ismax(j)
               nem=nem+1
               emspec(1,nem)=lbspec(j,i)
               emspec(2,nem)=fspec(j,i)
            enddo 
         else
            do i = 1,ismax(j)            
              if (lbspec(j,i).ge.lmin(j-1) .AND.
     >            lbspec(j,i).le.lmax(j-1)) then
                 do k=1,nem
                   if ( ABS(lbspec(j,i)-emspec(1,k)).lt.0.1) then
                      emspec(2,k)=emspec(2,k) + fspec(j,i)
                   endif
                 enddo
              elseif (lbspec(j,i).gt.lmax(j-1)) then
                 nem=nem+1
                 emspec(1,nem)=lbspec(j,i)
                 emspec(2,nem)=fspec(j,i)
              endif
            enddo
        endif 
      enddo 
      nem=nem+1
      emspec(1,nem) = 10000000.d0
      emspec(2,nem) = 0.d0
           
        
      if (nem.gt.2000) write(*,*) 'Pb of dimension for EM spectrum'
c
      return
      end
