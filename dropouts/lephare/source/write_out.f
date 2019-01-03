c
c     description : write string with output 
c
      subroutine WRITE_OUT(wpara,iwpara,spec,
     > zmin,chimin,imasmin,agemin,extilmin,extmin,zfmin,dmmin,
     > pdz,mag_abs,cont,zs,str_inp,zb,
     > chi,imasb,ageb,extb,zfb,dmb,mag_absb,nb,
     > ab,sab,kap,mabs,imagm,goodfilter,magm,
     > dzpdz,npdz,
     > z68i,z68s,z90i,z90s,z99i,z99s,
     > zml68i,zml68s,zml90i,zml90s,zml99i,zml99s,
     > zbay,zbayi,zbays,
     > mabsq,kapq,
     > nbused,nbul,
     > mod_distb,
     > parainf,parasup,paramed,
     > lumfirb,libfirb,modfirb,chi2_fir,nf_fir,dmfirb,
     > lirmed,lirinf,lirsup,absfir,kcfirb,magfirb,
     >  physpbest,
     >  ppbest,chipbest,reclpbest,
     >  ppmed,ppinf,ppsup,
     >  fluxphys,magphys0,kcorphys,
     >  zvmax,nzmax,
     > str_out,iwout)
c
      implicit none
c
      integer*4  chisize,nbf
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_zchi.decl'
      integer*4  i,j,k,pass,lnblnk,len
      integer*4  iwpara,iwout,jp
      integer*4  spec,imagm,npdz,nb(500)
      integer*4  imasmin(3),extilmin(3)
      integer*4  imasb(chisize),goodfilter(nbf)
      integer*4  nbused,nbul
c
      real*8     cont
      real*8     zmin(3),chimin(3),agemin(3),dmmin(3),extmin(3)
      real*8     zfmin(3),mag_abs(3),mod_distb(3)
      real*8     zb(chisize),chi(2,chisize),ageb(chisize),dmb(chisize)
      real*8     extb(chisize),zfb(chisize),mag_absb(chisize)
      real*8     pdz(500),zs,dzpdz(500)
      real*8     z68i,z68s,z90i,z90s,z99i,z99s
      real*8     zml68i,zml68s,zml90i,zml90s,zml99i,zml99s
      real*8     ab(nbf),sab(nbf),kap(nbf),mabs(nbf)
      real*8     kapq(nbf),mabsq(nbf),magm(nbf)
      real*8     zbay,zbayi,zbays
      integer*4  libfirb,modfirb,nf_fir
      real*8     chi2_fir,lumfirb,lirmed,lirinf,lirsup
      real*8     absfir(nbf),kcfirb(nbf),dmfirb,magfirb(nbf)       
      real*8     physpbest(50),fluxphys(nbf),magphys0(nbf),kcorphys(nbf)
      integer*4  reclpbest
      real*8     ppbest(50),ppmed(50),ppinf(50),ppsup(50),chipbest
      real*8     parainf(6),parasup(6),paramed(6) 
      integer*4  nzmax
      real*8     zvmax(500)

c
      character*4096  wpara(500),str_out(500),str,ich
      character*4096 str_inp,str2,paral
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c    initialisation
      do i = 1,500
         do k = 1,4096
            str_out(i)(k:k) = ' '
         enddo   
      enddo   
      do i = 1,4096
         str(i:i) = ' '
      enddo   
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c    printing each parameter
      k = 1
      do i = 1, iwpara
         pass = 0 
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c write output with specific format 
cccccccccc OBJECT  cccccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'IDENT' ) then 
            write(str,'(i10)') spec            
            pass = pass + 1
         endif   
ccccccccc  BEST FIT cccccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_BEST' ) then
            write(str,'(f9.4)') zmin(1)
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_BEST68_LOW' ) then
            write(str,'(f9.4)') z68i
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_BEST68_HIGH' ) then
            write(str,'(f9.4)') z68s
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_BEST90_LOW' ) then
            write(str,'(f9.4)') z90i
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_BEST90_HIGH' ) then
            write(str,'(f9.4)') z90s
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_BEST99_LOW' ) then
            write(str,'(f9.4)') z99i
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_BEST99_HIGH' ) then
            write(str,'(f9.4)') z99s
            pass = pass + 1
         endif   
c
c
c         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML68_LOW' ) then
c            write(str,'(f9.4)') zml68i
c            pass = pass + 1
c         endif   
c
c         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML68_HIGH' ) then
c            write(str,'(f9.4)') zml68s
c            pass = pass + 1
c         endif   
c
c         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML90_LOW' ) then
c            write(str,'(f9.4)') zml90i
c            pass = pass + 1
c         endif   
c
c         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML90_HIGH' ) then
c            write(str,'(f9.4)') zml90s
c            pass = pass + 1
c         endif   
c
c         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML99_LOW' ) then
c            write(str,'(f9.4)') zml99i
c            pass = pass + 1
c         endif   
c
c         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML99_HIGH' ) then
c            write(str,'(f9.4)') zml99s
c            pass = pass + 1
c         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML' ) then
            write(str,'(f9.4)') zbay
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML68_LOW' ) then
            write(str,'(f9.4)') zbayi
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML68_HIGH' ) then
            write(str,'(f9.4)') zbays
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'CHI_BEST' ) then
            write(str,'(E12.6)') chimin(1) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MOD_BEST' ) then
            write(str,'(i6)') imasmin(1) 
            pass = pass + 1
         endif   
c
c         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'AGE_BEST' ) then
c            write(str,'(E12.6)') agemin(1) 
c            pass = pass + 1
c         endif
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'EXTLAW_BEST' ) then
            write(str,'(i6)') extilmin(1) 
            pass = pass + 1
         endif   
   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'EBV_BEST' ) then
            write(str,'(f8.3)') extmin(1) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ZF_BEST' ) then
            write(str,'(f9.4)') zfmin(1) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'SCALE_BEST' ) then
            write(str,'(E12.6)') dmmin(1) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PDZ_BEST' ) then
            write(str,'(f8.3)') pdz(1) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_ABS_BEST' ) then
            write(str,'(f9.3)') mag_abs(1) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIST_MOD_BEST' ) then
            write(str,'(E12.6)') mod_distb(1) 
            pass = pass + 1
         endif   

ccccccccc   INPUT INFO ccccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'CONTEXT' ) then
            write(str,'(g20.7)') cont 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ZSPEC' ) then
            write(str,'(E12.6)') zs 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'STRING_INPUT' ) then
            str = str_inp
            pass = pass + 1
         endif   
cccccccccc  BEST   QSO     ccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_QSO' ) then
            write(str,'(f9.4)') zmin(2)
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'CHI_QSO' ) then
            write(str,'(E12.6)') chimin(2) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MOD_QSO' ) then
            write(str,'(i6)') imasmin(2) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_ABS_QSO' ) then
            write(str,'(f9.3)') mag_abs(2) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIST_MOD_QSO' ) then
            write(str,'(E12.6)') mod_distb(2) 
            pass = pass + 1
         endif   
ccccccccc  BEST   STAR   cccccccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MOD_STAR' ) then
            write(str,'(i6)') imasmin(3) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'CHI_STAR' ) then
            write(str,'(E12.6)') chimin(3) 
            pass = pass + 1
         endif   
ccccccccc    SECONDARY GALAXY PEAK   cccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_SEC' ) then
            write(str,'(f9.4)') zb(nb(2))
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'CHI_SEC' ) then
            write(str,'(E12.6)') chi(2,nb(2)) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MOD_SEC' ) then
            write(str,'(i6)') imasb(nb(2)) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'AGE_SEC' ) then
            write(str,'(E12.6)') ageb(nb(2)) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'EBV_SEC' ) then
            write(str,'(f8.3)') extb(nb(2)) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ZF_SEC' ) then
            write(str,'(f9.4)') zfb(nb(2)) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'SCALE_SEC' ) then
            write(str,'(E12.6)') dmb(nb(2)) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PDZ_SEC' ) then
            write(str,'(f8.3)') pdz(2) 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_ABS_SEC' ) then
            write(str,'(f9.3)') mag_absb(nb(2)) 
            pass = pass + 1
         endif            
cccccccccccc   OBSERVED  AB  MAGNITUDES   ccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_OBS()' ) then
c           write(str,'(50(e12.6,1x))') (ab(j),j=1,imagm) 
           write(str,'(100(f9.3,1x))') (ab(j),j=1,imagm) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_OBS()' ) then
c           write(str,'(50(e12.6,1x))') (sab(j),j=1,imagm) 
           write(str,'(100(f9.3,1x))') (sab(j),j=1,imagm) 
           pass = pass + 1
         endif  
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_MOD()' ) then
           write(str,'(100(f9.3,1x))') (magm(j),j=1,imagm) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PHYS_MAG_MOD()' ) then
           write(str,'(100(f9.3,1x))') (fluxphys(j),j=1,imagm) 
           pass = pass + 1
         endif  
c
cccccccccccc  K-CORRECTIONS in ALL BANDS  cccccccccccc 
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'K_COR()' ) then
           write(str,'(100(f9.3,1x))') (kap(j),j=1,imagm) 
           pass = pass + 1
         endif 
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PHYS_K_COR()' ) then
           write(str,'(100(f9.3,1x))') (kcorphys(j),j=1,imagm) 
           pass = pass + 1
         endif 
c
cccccccccccc  ABSOLUTE MAGNITUDES in ALL BANDS ccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_ABS()' ) then
           write(str,'(100(f9.3,1x))') (mabs(j),j=1,imagm) 
           pass = pass + 1
         endif  
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PHYS_MAG_ABS()' ) then
           write(str,'(100(f9.3,1x))') (magphys0(j),j=1,imagm) 
           pass = pass + 1
         endif  

cccccccccccc  ADOPTED FILTER FOR MAG_ABS in ALL BANDS ccccccc
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MABS_FILT()' ) then
            write(str,'(100(i6,1x))') (goodfilter(j),j=1,imagm) 
            pass = pass + 1 
         endif   
         
cccccccccccc  K-CORRECTIONS in ALL BANDS for QSO cccccccccccc 
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'K_COR_QSO()' ) then
           write(str,'(100(f9.3,1x))') (kapq(j),j=1,imagm) 
           pass = pass + 1
         endif  
cccccccccccc  ABSOLUTE MAGNITUDES in ALL BANDS for QSO ccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_ABS_QSO()' ) then
           write(str,'(100(f9.3,1x))') (mabsq(j),j=1,imagm) 
           pass = pass + 1
         endif  
cccccccccccc  REDSHIFT PROBABILITY IN Z-BINS  cccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PDZ()' ) then
            if (npdz.gt.1) then 
              write(str,'(500(f9.4,1x))') (dzpdz(j),j=1,npdz/2) 
              pass = pass + 1
            endif
         endif  

cccccccccccc  NUMBER OF BANDS USED   cccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'NBAND_USED' ) then
            write(str,'(i4)') nbused 
            pass = pass + 1
         endif   
cccccccccccc  NUMBER OF BANDS USED WITH UPPER LIMITS  cccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'NBAND_ULIM' ) then
            write(str,'(i4)') nbul 
            pass = pass + 1
         endif   
cccccccccccc  FIR OUPUT  
c       integer*4  libfirb,modfirb,nf_fir
c       real*8     chi2_fir,lumfirb,lirmed,lirinf,lirsup       
        if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LUM_TIR_BEST' ) then
            write(str,'(f9.4)')  lumfirb
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LUM_TIR_MED'  ) then
            write(str,'(f9.4)')  lirmed
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LUM_TIR_INF'  ) then
            write(str,'(f9.4)')  lirinf
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LUM_TIR_SUP'  ) then
            write(str,'(f9.4)')  lirsup
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'CHI2_FIR' ) then
            write(str,'(E13.6)')  chi2_fir
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'FSCALE_FIR' ) then
            write(str,'(E13.6)')  dmfirb
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'NBAND_FIR' ) then
            write(str,'(i6)') nf_fir 
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LIB_FIR' ) then
            write(str,'(i6)') libfirb 
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MOD_FIR' ) then
            write(str,'(i6)') modfirb
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_MOD_FIR()' ) then
           write(str,'(100(f9.3,1x))') (magfirb(j),j=1,imagm) 
           pass = pass + 1
         endif  
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_ABS_FIR()' ) then
           write(str,'(100(f9.3,1x))') (absfir(j),j=1,imagm) 
           pass = pass + 1
         endif  
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'K_COR_FIR()' ) then
           write(str,'(100(f9.3,1x))') (kcfirb(j),j=1,imagm) 
           pass = pass + 1
         endif 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Physical parameters from galaxy library
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MASS_BEST' ) then
               write(str,'(E13.6)')  physpbest(6)
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'SFR_BEST' ) then
               write(str,'(E13.6)')  physpbest(7)
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'SSFR_BEST' ) then
               write(str,'(E13.6)')  physpbest(7)-physpbest(6)
            pass = pass + 1
         endif   

         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LDUST_BEST' ) then
               write(str,'(E13.6)')  physpbest(5)
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'AGE_BEST' ) then
               write(str,'(E13.6)')  physpbest(1)
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LUM_NUV_BEST' ) then
               write(str,'(E13.6)')  physpbest(2)
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LUM_R_BEST' ) then
               write(str,'(E13.6)')  physpbest(3)
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LUM_K_BEST' ) then
               write(str,'(E13.6)')  physpbest(4)
            pass = pass + 1
         endif   
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_MAX()' ) then
               write(str,'(500(E13.6,2x))')  (zvmax(j),j=1,nzmax)
            pass = pass + 1
         endif   

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Phys parameters from ZPHOTLIB
         do jp = 1, 6
            if (jp==1)  paral = "AGE"
            if (jp==2)  paral = "LDUST"
            if (jp==3)  paral = "LDUST2"
            if (jp==4)  paral = "MASS"
            if (jp==5)  paral = "SFR"
            if (jp==6)  paral = "SSFR"
c            
            str2= paral(1:lnblnk(paral)) //"_MED"            
            if (wpara(i)(1:lnblnk(wpara(i))) .eq. 
     >                       str2(1:lnblnk(str2))) then
                write(str,'(E13.6)')  paramed(jp)
               pass = pass + 1
            endif   
c            
            str2= paral(1:lnblnk(paral)) //"_INF"            
            if (wpara(i)(1:lnblnk(wpara(i))) .eq. 
     >                       str2(1:lnblnk(str2))) then
                write(str,'(E13.6)')  parainf(jp)
               pass = pass + 1
            endif   
c            
            str2= paral(1:lnblnk(paral)) //"_SUP"            
            if (wpara(i)(1:lnblnk(wpara(i))) .eq. 
     >                       str2(1:lnblnk(str2))) then
                write(str,'(E13.6)')  parasup(jp)
               pass = pass + 1
            endif   
         enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      PHYSICAL PARAMETERS FROM STOCHASTIC LIBRARY (LIB_PHYS)
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PHYS_MOD_BEST' ) then
            write(str,'(i12)') reclpbest-1 
            pass = pass + 1 
         endif           
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PHYS_CHI2_BEST' ) then
           write(str,'(E12.6)') chipbest 
           pass = pass + 1
         endif  
c

         if (wpara(i)(1:9) .eq. 'PHYS_PARA')then
          do jp = 1, 28 
           if (jp.le.9)  write(ich,'(i1.1)') jp  
           if (jp.ge.10) write(ich,'(i2.2)') jp 

           str2= "PHYS_PARA"// ich(1:lnblnk(ich)) //"_BEST"
           if(wpara(i)(1:lnblnk(wpara(i))).eq.str2(1:lnblnk(str2)))then
             write(str,'(E13.6)') ppbest(jp) 
             pass = pass + 1
           endif  
c
           str2= "PHYS_PARA"// ich(1:lnblnk(ich)) //"_MED"
           if(wpara(i)(1:lnblnk(wpara(i))).eq.str2(1:lnblnk(str2)))then
              write(str,'(E13.6)') ppmed(jp)
              pass = pass + 1
           endif  

           str2= "PHYS_PARA"// ich(1:lnblnk(ich)) //"_INF"
           if(wpara(i)(1:lnblnk(wpara(i))).eq.str2(1:lnblnk(str2)))then
              write(str,'(E13.6)') ppinf(jp)
              pass = pass + 1
           endif  

           str2= "PHYS_PARA"// ich(1:lnblnk(ich)) //"_SUP"
           if(wpara(i)(1:lnblnk(wpara(i))).eq.str2(1:lnblnk(str2)))then
              write(str,'(E13.6)') ppsup(jp)
              pass = pass + 1
           endif  
          enddo
         endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (pass .eq. 1) then 
c test length for the output 
           len = lnblnk(str_out(k))+lnblnk(str) +2 
           if (len.ge.4096) k = k + 1
c print output in str_out 
           str_out(k) = str_out(k)(1:lnblnk(str_out(k)))
     >     // ' ' // str(1:lnblnk(str))
c
c           write(*,*) i,k,str_out(k)(1:lnblnk(str_out(k)))
         endif  
      enddo
      iwout = k    
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
