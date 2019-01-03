c
c     description : write string with output 
c
      subroutine WRITE_SCOUT(wpara,iwpara,spec,
     > zmin,chimin,imasmin,agemin,extmin,zfmin,dmmin,pdz,mag_abs,
     > cont,zs,str_inp,zb,
     > chi,imasb,ageb,extb,zfb,dmb,mag_absb,nb,
     > ab,sab,kap,mabs,imagm,goodfilter,magm,
     > dzpdz,npdz,
     > z68i,z68s,z90i,z90s,z99i,z99s,
     > zml68i,zml68s,zml90i,zml90s,zml99i,zml99s,
     > zbay,zbayi,zbays,
     > mabsq,kapq,
     > nbused,nbul,
     > str_out,iwout,
     > parasc
     > )
c
      implicit none
c
      integer*4  chisize,nbf
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_zchi.decl'
      integer*4  i,j,k,pass,lnblnk,len
      integer*4  iwpara,iwout
      integer*4  spec,cont,imagm,npdz,nb(100)
      integer*4  imasmin(3),imasb(chisize),goodfilter(nbf)
      integer*4  nbused,nbul
c
      real*8     zmin(3),chimin(3),agemin(3),dmmin(3),extmin(3)
      real*8     zfmin(3),mag_abs(3)
      real*8     zb(chisize),chi(2,chisize),ageb(chisize),dmb(chisize)
      real*8     extb(chisize),zfb(chisize),mag_absb(chisize)
c      real*8     mod_distb(3),dist_mod(chisize)
      real*8     pdz(100),zs,dzpdz(50)
      real*8     z68i,z68s,z90i,z90s,z99i,z99s
      real*8     zml68i,zml68s,zml90i,zml90s,zml99i,zml99s
      real*8     ab(nbf),sab(nbf),kap(nbf),mabs(nbf)
      real*8     kapq(nbf),mabsq(nbf),magm(nbf)
      real*8     zbay,zbayi,zbays
      real*8     parasc(10,3)
c
      character*512  wpara(100),str_out(100),str
      character*1024 str_inp
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c    initialisation
      do i = 1,100
         do k = 1,512
            str_out(i)(k:k) = ' '
         enddo   
      enddo   
      do i = 1,512
         str(i:i) = ' '
      enddo   
c
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
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML68_LOW' ) then
            write(str,'(f9.4)') zml68i
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML68_HIGH' ) then
            write(str,'(f9.4)') zml68s
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML90_LOW' ) then
            write(str,'(f9.4)') zml90i
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML90_HIGH' ) then
            write(str,'(f9.4)') zml90s
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML99_LOW' ) then
            write(str,'(f9.4)') zml99i
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_ML99_HIGH' ) then
            write(str,'(f9.4)') zml99s
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_BAYSIAN' ) then
            write(str,'(f9.4)') zbay
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_BAYS68_LOW' ) then
            write(str,'(f9.4)') zbayi
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_BAYS68_HIGH' ) then
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
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'AGE_BEST' ) then
            write(str,'(E12.6)') agemin(1) 
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
ccccccccc   INPUT INFO ccccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'CONTEXT' ) then
            write(str,'(i9)') cont 
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ZSPEC' ) then
            write(str,'(f9.4)') zs 
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
           write(str,'(100(f9.3,1x))') (ab(j),j=1,imagm) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_OBS()' ) then
           write(str,'(100(f9.3,1x))') (sab(j),j=1,imagm) 
           pass = pass + 1
         endif  
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_MOD()' ) then
           write(str,'(100(f9.3,1x))') (magm(j),j=1,imagm) 
           pass = pass + 1
         endif  
c
cccccccccccc  K-CORRECTIONS in ALL BANDS  cccccccccccc 
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'K_COR()' ) then
           write(str,'(100(f9.3,1x))') (kap(j),j=1,imagm) 
           pass = pass + 1
         endif  
cccccccccccc  ABSOLUTE MAGNITUDES in ALL BANDS ccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_ABS()' ) then
           write(str,'(100(f9.3,1x))') (mabs(j),j=1,imagm) 
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
              write(str,'(100(f9.4,1x))') (dzpdz(j),j=1,npdz/2) 
              pass = pass + 1
            endif
         endif  

cccccccccccc  PARAMETERS from SC LIBRARY    ccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PARASC_MED()') then
           write(str,'(100(f9.3,1x))') (parasc(k,1),k=1,7) 
           pass = pass + 1
         endif  
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PARASC68_LOW()') then
           write(str,'(100(f9.3,1x))') (parasc(k,2),k=1,7) 
           pass = pass + 1
         endif  
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PARASC68_HIGH()') then
           write(str,'(100(f9.3,1x))') (parasc(k,3),k=1,7) 
           pass = pass + 1
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (pass .eq. 1) then 
c test length for the output 
           len = lnblnk(str_out(k))+lnblnk(str) +2 
           if (len.ge.512) k = k + 1
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
