c
c     description : write string with output for simulation 
c
      subroutine WRITE_SIM(wpara,iwpara,sb,nobj,
     >  mags,zref,lf,kciref,ageref,ebvref,zforref,
     >  mabs,rnum,ecor,phicor,alpcor,
     >  nfilt,magrefo,mabsz,re0,re0as,mabs_bv,thetamag,
     >  magref,errmag,flran,errflran,magsb,errmagsb,dl,
     > str_out,iwout)
c
      implicit none
c
      integer*4  nbf
      INCLUDE 'dim_filt.decl'
      integer*4  i,j,k,pass,lnblnk,len
      integer*4  iwpara,iwout
      integer*4  nobj,lf,rnum,nfilt
c      integer*4  
c
      real*8    mags,zref,kciref,mabs,dl
      real*8    ecor,phicor,alpcor 
      real*8    re0,re0as,mabs_bv,ageref,ebvref,zforref 
      real*8    magrefo(nbf),magref(nbf),errmag(nbf)
      real*8    flran(nbf),errflran(nbf)
      real*8    mabsz(nbf) 
      real*8    magsb(10,nbf),errmagsb(10,nbf) 
      real*8    thetamag(10,nbf) 
c
      character*512 wpara(100),str_out(100),str,sb
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
            write(str,'(i10)') nobj            
            pass = pass + 1
         endif   
ccccccccc  REFERENCE INFO cccccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'Z_REF' ) then
            write(str,'(f8.3)') zref
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MODEL_REF' ) then
            write(str,'(i8)') rnum
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'AGE_REF' ) then
            write(str,'(E12.6)') ageref
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'EBV_REF' ) then
            write(str,'(f8.3)') ebvref
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ZFOR_REF' ) then
            write(str,'(f8.3)') zforref
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LF_REF' ) then
            write(str,'(i8)') lf
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_REF' ) then
            write(str,'(f8.3)') mags
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MABS_REF' ) then
            write(str,'(f8.3)') mabs
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'KCOR_REF' ) then
            write(str,'(f9.3)') kciref
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'LUM_EVOL' ) then
            write(str,'(f9.3)') ecor
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'PHI_EVOL' ) then
            write(str,'(E12.6)') phicor
            pass = pass + 1
         endif   
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ALPHA_EVOL' ) then
            write(str,'(f8.3)') alpcor
            pass = pass + 1
         endif   
c
          if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIST_LUM' ) then
            write(str,'(E12.6)') dl 
            pass = pass + 1
          endif
c
         if (sb(1:1).eq.'y' .or. sb(1:1).eq.'Y') then 

          if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'RE_REF' ) then
            write(str,'(E12.6)') re0
            pass = pass + 1
          endif   
c
          if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'RE_ARCSEC' ) then
            write(str,'(E12.6)') re0as
            pass = pass + 1
          endif
c
          if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'BABS_REF' ) then
            write(str,'(f8.3)') mabs_bv
            pass = pass + 1
          endif   
         endif
cccccccccccc   MAGNITUDES WITH  NO SB MODE    ccccccccccc
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_REF()' ) then
           write(str,'(100(f9.3,1x))') (magrefo(j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_RAN()' ) then
           write(str,'(100(f9.3,1x))') (magref(j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_RAN()' ) then
           write(str,'(100(f9.3,1x))') (errmag(j),j=1,nfilt) 
           pass = pass + 1
         endif  
cccccccccccc   FLUXES WITH  NO SB MODE    ccccccccccc
c         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'FLUX_REF()' ) then
c           write(str,'(100(E12.6,1x))') (magrefo(j),j=1,nfilt) 
c           pass = pass + 1
c         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'FLUX_RAN()' ) then
           write(str,'(100(E12.6,1x))') (flran(j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_FLUX_RAN()' ) then
           write(str,'(100(E12.6,1x))') (errflran(j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
cccccccccccc   MAGNITUDES WITH  SB MODE    ccccccccccc
         if (sb(1:1).eq.'y' .or. sb(1:1).eq.'Y') then 
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_ISO()' ) then
           write(str,'(100(f9.3,1x))') (magsb(1,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_ISO()' ) then
           write(str,'(100(f9.3,1x))') (errmagsb(1,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIAM_ISO()' ) then
           write(str,'(100(f9.3,1x))') (thetamag(1,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_PSEUDO()' ) then
           write(str,'(100(f9.3,1x))') (magsb(2,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_PSEUDO()' ) then
           write(str,'(100(f9.3,1x))') (errmagsb(2,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIAM_PSEUDO()' ) then
           write(str,'(100(f9.3,1x))') (thetamag(2,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_AUTO()' ) then
           write(str,'(100(f9.3,1x))') (magsb(3,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_AUTO()' ) then
           write(str,'(100(f9.3,1x))') (errmagsb(3,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIAM_AUTO()' ) then
           write(str,'(100(f9.3,1x))') (thetamag(3,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_PETRO()' ) then
           write(str,'(100(f9.3,1x))') (magsb(4,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_PETRO()' ) then
           write(str,'(100(f9.3,1x))') (errmagsb(4,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIAM_PETRO()' ) then
           write(str,'(100(f9.3,1x))') (thetamag(4,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_APER1()' ) then
           write(str,'(100(f9.3,1x))') (magsb(5,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_APER1()' ) then
           write(str,'(100(f9.3,1x))') (errmagsb(5,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIAM_APER1()' ) then
           write(str,'(100(f9.3,1x))') (thetamag(5,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_APER2()' ) then
           write(str,'(100(f9.3,1x))') (magsb(6,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_APER2()' ) then
           write(str,'(100(f9.3,1x))') (errmagsb(6,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIAM_APER2()' ) then
           write(str,'(100(f9.3,1x))') (thetamag(6,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_APER3()' ) then
           write(str,'(100(f9.3,1x))') (magsb(7,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_APER3()' ) then
           write(str,'(100(f9.3,1x))') (errmagsb(7,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIAM_APER3()' ) then
           write(str,'(100(f9.3,1x))') (thetamag(7,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_APER4()' ) then
           write(str,'(100(f9.3,1x))') (magsb(8,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_APER4()' ) then
           write(str,'(100(f9.3,1x))') (errmagsb(8,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIAM_APER4()' ) then
           write(str,'(100(f9.3,1x))') (thetamag(8,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'MAG_APER5()' ) then
           write(str,'(100(f9.3,1x))') (magsb(9,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'ERR_MAG_APER5()' ) then
           write(str,'(100(f9.3,1x))') (errmagsb(9,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
         if (wpara(i)(1:lnblnk(wpara(i))) .eq. 'DIAM_APER5()' ) then
           write(str,'(100(f9.3,1x))') (thetamag(9,j),j=1,nfilt) 
           pass = pass + 1
         endif  
c
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
