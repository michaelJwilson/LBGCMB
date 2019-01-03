c   last modif : 10/02/03
      program mag_zform
c     1 option : mag_zform  
c                        -c zphot.para    : config file
c     
c     Construction des differents modeles d'evolution pour
c     la mesure des z-photometriques fonction de Omo,Oml,Ho des
c     bandes photometriques et de l'extinction
c
c     Note :  Only one zformation per model 
c             Mag , kcor  are interpolated to fit with z,dz 
c             the value for extinction is lost and replaced by zform 
c             only k+e correction writen in binary output file 
c              ---> The library can be used for zphot !  
c             The output format is : 
c             model,extinction,z,age,zform,mag,kcor, k+e-cor, kappa-cor
      implicit none
      integer*4  nbf,maxsize
      integer*4  wmax,amax
      INCLUDE 'dim_wave.decl'      
c      parameter  (wmax=8000, amax=500)
      parameter  (amax=500)
      parameter  (maxsize=110000)
      INCLUDE 'dim_filt.decl'
      integer*4  jmax(nbf),pass    
      integer*4  imag,iass,imas,i,j,k,n,ie,if,iw,nebv,nzf
      integer*4  m,iassm,test,agemax
      integer*4  recmax,nrec(maxsize),reclmax,irec
      integer*4  nmod,nr,jtot,orec,orecmax,nlveg,lnblnk
      integer*4  paravi(500),modext(2),iext(10),iopa(81),wnmax
      integer*4  iextinc
c
      real*8     h0,om0,l0,dz,zmax,dzup,val,dummy,funz
      real*8     timy,eta,zoss,distmod0,distmodz
      real*8     wave(wmax),exti(10,2,wmax),extinc(2,wmax)
      real*8     vec(wmax),zf(500)
      real*8     age(amax),fl(amax,wmax),flux(wmax),fext(wmax)
      real*8     paravr(500),ebv(500),mag(nbf),flux0(wmax)
      real*8     ecor(nbf),kcor(nbf),kecor(nbf),kapcor(nbf)
      real*8     mag0(nbf),mag1(nbf),mag2(nbf)
      real*8     lbveg(wmax),flveg(wmax)
      real*8     lambf(nbf,maxsize), repf(nbf,maxsize)
      real*8     opal(81,wmax),opat(81,wmax)
      real*8     wnew(wmax),fnew(wmax)
      real*8     zp(nbf),flmoy(nbf),flwidth(nbf),abcor(nbf)
      real*8     zform,tform,toss,zeta,z0
      real*8     fcorr(nbf)
c      
      integer*4  iz,izl,ma
      real*8     fluxs(wmax),fluxi(wmax),zs,zi
      real*8     magi(nbf),mags(nbf),mag2i(nbf),mag2s(nbf)
      real*8     mag0i(nbf),mag0s(nbf)
      real*8     distmodzi,distmodzs
      real*8     ecori,kcori,kapcori
      real*8     ecors,kcors,kapcors
      real*8     interp
c
      character*4096  name(nbf),magtyp,param,vega
      character*4096  paravc(500),str,zpwork
      character*4096  valc,fileop(500),sedtyp
      character*4096  file,filters,lib,colib,zpdir,extlaw,config,outasc
      character*4096 cr_date,fdate
      external   timy,funz,zeta,interp
c
      common /univ/ h0,om0,l0
      common /pas/ dz
      common /input/ iass, imas
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc   PARAMETERS                         cccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  environmental variable 
      call getenv('LEPHAREDIR',zpdir)
      test=lnblnk(zpdir)
      if (test .eq. 0) then
        write(6,*) 'WARNING :  variable LEPHAREDIR not defined'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(6,*) 'WARNING :  variable LEPHAREWORK not defined'
        stop
      endif
cccccccccccccccccccccccccccccccccccccccccccccccc
c help on line
      param='mag_zform'
      call get_help(test)
      if (test .eq. 1) call help(param)
cccccccccccccccccccccccccccccccccccccccccccccccc
c Initialisation of Input parameter from config file
      param='-c'
      call get_conf(param,config,test)
      if (test.ne.1) call err_option(param,1)
      call get_path(config)
c      param='-t'
c      call get_conf(param,sedtyp,test)
c      call getc_option(param,config,1,paravc,test)
c      if (test.eq.1) sedtyp = paravc(1)(1:lnblnk(paravc(1)))
c      if (test.ne.1) call err_option(param,1)
      sedtyp='GAL'
c  read option
c
      param='-COSMOLOGY'
      call getf_option(param,config,3,paravr,test)
      if (test.eq.3)   h0=paravr(1)
      if (test.eq.3)   om0=paravr(2)
      if (test.eq.3)   l0=paravr(3)
      if (test.ne.3)  call err_option(param,1)
c
      param='-FILTER_FILE'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  filters=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
c
      param='-MAGTYPE'
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1) magtyp=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
c
      param='-EXTINC_LAW'
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1) extlaw=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  then
        call err_option(param,2)
        extlaw='NONE'
      endif  
c
      if (extlaw(1:4).eq.'NONE' .or. extlaw(1:4).eq.'none') then
         nebv=1
         ebv(1)=0.
         modext(1) = 0
         modext(2) = 0
      else   
         param='-EB_V'
         call getf_option(param,config,100,paravr,test)
         nebv=test
         if (test.ge.1) then
            do i = 1,nebv
              ebv(i) =   paravr(i)
            enddo
         else
           nebv=1
           ebv(1)=0
         endif   
         if (test.eq.0)  call err_option(param,2)
         param='-MOD_EXTINC'
         call geti_option(param,config,2,paravi,test)
         if (test.eq.2)  modext(1) = paravi(1)
         if (test.eq.2)  modext(2) = paravi(2)
         if (test.ne.2)  modext(1) = 0
         if (test.ne.2)  modext(2) = 0 
         if (test.ne.2)  call err_option(param,2)
      endif   
c
      param='-Z_FORM'
      call getf_option(param,config,100,paravr,test)
      nzf=test
      if (test.ge.1) then
         do i = 1,nzf
           zf(i) =   paravr(i)
         enddo
      else
         nzf=1
         zf(1)=0
      endif   
      if (test.eq.0)  call err_option(param,1)      
c
      param='-Z_STEP'
      call getf_option(param,config,3,paravr,test)
      if (test.eq.3)  dz = paravr(1)
      if (test.eq.3)  zmax = paravr(2)
      if (test.eq.3)  dzup = paravr(3)
      if (test.ne.3)  call err_option(param,1)
c
      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
        param='-GAL_LIB_IN'
      elseif (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q") then
        param='-QSO_LIB_IN'
      endif 
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1)  lib=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
c
      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
         param='-GAL_LIB_OUT'
      elseif (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q") then
         param='-QSO_LIB_OUT'
      endif
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1)  colib=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
c
      param='-LIB_ASCII'
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1)  outasc=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,2)
      if (test.ne.1)  outasc='NO'
c
c  writing options    
       write(6,'(A)')"############################################"
       write(6,'(A)')"# It s computing the SYNTHETIC MAGNITUDES  #"
       write(6,'(A)')"# For Gal/QSO libraries with these OPTIONS #"
       file = zpwork(1:lnblnk(zpwork)) //'/filt/'
     >      // filters(1:lnblnk(filters))
       write(6,'(2A)')"# FILTER_FILE : ",file(1:lnblnk(file))  
       write(6,'(2A)')"# MAGTYPE     : ",magtyp(1:lnblnk(magtyp))
       if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
           file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >      // lib(1:lnblnk(lib)) //'(.doc & .bin)'
         write(6,'(2A)') "# GAL_LIB_IN  : ",file(1:lnblnk(file)) 
         file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >      // colib(1:lnblnk(colib)) //'(.doc & .bin)'
         write(6,'(2A)') "# GAL_LIB_OUT : ",file(1:lnblnk(file))     
         write(6,'(A,2(I6,1x))') "# MOD_EXTINC  : ",modext(1),modext(2)
       elseif (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q") then
           file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >      // lib(1:lnblnk(lib)) //'(.doc & .bin)'
         write(6,'(2A)') "# QSO_LIB_IN  : ",file(1:lnblnk(file)) 
          file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >      // colib(1:lnblnk(colib)) //'(.doc & .bin)'
         write(6,'(2A)') "# QSO_LIB_OUT : ",file(1:lnblnk(file))     
       endif  
       write(6,'(A,50(f5.2,1x))') "# Z_FORM     : ",(zf(i),i=1,nzf) 
       write(6,'(A,3(f7.3,1x))') "# Z_STEP      : ",dz,zmax,dzup    
       write(6,'(A,3(f7.3,1x))') "# COSMOLOGY   : ",h0,om0,l0
       write(6,'(2A)') "# EXTINC_LAW  : ",extlaw(1:lnblnk(extlaw))   
       write(6,'(A,10(f6.3,1x))') "# EB_V        :",(ebv(i),i=1,nebv)
       write(6,'(2A)') "# LIB_ASCII   : ",outasc(1:lnblnk(outasc)) 
       write(6,'(A)')"############################################"
c
ccccccccccccccccccccccccccccccccccc
c      READING FILES  
ccccccccccccccccccccccccccccccc
c  read doc file from GAL_LIB_IN 
      file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     > // lib(1:lnblnk(lib)) // '.doc'
      param='NUMBER_SED'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') imas
      do i = 1, imas
        write(str,'(i6.6)') i 
         param = 'MOD_' // str(1:lnblnk(str))
         call read_para2(param,file,paravc,test)
         if (test.ne.0) then
            read(paravc(3),'(i8)') nrec(i)
         else
            call err_option(param,1)
         endif   
      enddo   
      param='NUMBER_ROWS'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') recmax
      param='RECORD_LENGTH'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') reclmax      
c
c  check if nzf = nb models
      if (nzf .ne. imas) then
         write(6,*) ' Nb Zform differs fron Nb models --> STOP'
         stop
      endif   
c  read Vega spectrum
      vega = zpdir(1:lnblnk(zpdir)) // '/vega/VegaLCB.sed'
      open(1, file=vega,status='old')
      read(1,*)   
      i=0
      do while (.true.)
	i=i+1
        read(1,*,end=10)  lbveg(i), flveg(i)
      enddo
 10   nlveg=i-1  
      close(1)
c  read the filters from  file
      file = zpwork(1:lnblnk(zpwork)) //'/filt/'
     >        // filters(1:lnblnk(filters))
      open(1,file=file)
      read(1,*) valc,imag
      do i =  1, imag
         read(1,*) valc,jmax(i),name(i)
         do j = 1,jmax(i)
	    read(1,*)  lambf(i,j), repf(i,j)
 	 enddo
      enddo
      close(1) 
c  getting AB corrections for info 
      call zeropoint(filters,zp,abcor,flmoy,flwidth,fcorr,imag)
      iext(1)=0
      if (extlaw(1:4).ne.'NONE' .and. extlaw(1:4).ne.'none') then 
c  read extinction law 
        file = zpdir(1:lnblnk(zpdir)) //'/ext/' 
     > // extlaw(1:lnblnk(extlaw))
        open(1,file=file)
        i=0
        do while (.true.)
          i=i+1
          read(1,*,end=11) exti(1,1,i),exti(1,2,i)
        enddo
 11     iext(1)=i-1
        close(1)
      endif  
c  reading file with  extragalacitic opacity
      file = zpdir(1:lnblnk(zpdir)) // '/opa/OPACITY.dat'
      open(1,file=file)
      do i = 1,81
          read(1,*) val, fileop(i)
      enddo
      close(1)
      do i = 1, 81
          file = zpdir(1:lnblnk(zpdir)) // '/opa/' 
     > // fileop(i)
          open(1,file=file)
          k = 0
          do while (.true.)
             k = k + 1
             read(1,*,end=12) opal(i,k), opat(i,k)
          enddo  
 12       iopa(i) = k - 1
          close(1)
      enddo    
      write(6,650) imag,(name(j),j=1,imag)
 650  format('* ',I4,' filters selected : ',100(A10,1x))
      write(6,900) imas
 900  format('* number of  models used =',i3)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c    OPENING INPUT and OUTPUT FILES           ccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Open input library with models ...
      file = zpwork(1:lnblnk(zpwork)) // '/lib_bin/'
     >       // lib(1:lnblnk(lib)) // '.bin'
      open(1,file=file,status='unknown',access='direct',recl=reclmax)
c  Open output mag lib file
c      orecmax = 8*(5+imag*2)   ! if kcorrection added 
      orecmax = 8*(6+imag*2)   ! if kcorrection added 
c      orecmax = 4*(5+imag)
      file = zpwork(1:lnblnk(zpwork)) // '/lib_mag/'
     >    // colib(1:lnblnk(colib)) // '.bin'
      open(2,file=file,form='unformatted',status='unknown',
     > access='direct',recl=orecmax)
c  Open output doc file
      file = zpwork(1:lnblnk(zpwork)) // '/lib_mag/'
     >    //colib(1:lnblnk(colib)) // '.doc'
      open(4,file=file)
c Open output mag ASCII file
      if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y") then 
        file = colib(1:lnblnk(colib)) // '.dat'
        open(5,file=file)
        write(5,'(A,1x,I6,1x,A,1x,20(A,2x))') "#",imag,magtyp(1:5),
     > (name(j)(1:lnblnk(name(j))),j=1,imag)
      endif
cccccccccccccccccccccccccccccccccccccccccccccc
c   Computing Magnitudes        cccccccccccccc   
cccccccccccccccccccccccccccccccccccccccccccccc
      irec = 0
      orec = 0
      pass = 0
      agemax=1
      tform=0
      do i = 1, imas                      ! loop on  SEDs 
c Read flux and lambda for each SED and number of age (n)
         n = 0
         do j = 1, nrec(i)+1
            irec = irec + 1
            read(1,rec=irec) nr,nmod,dummy,iw,(vec(k),k=1,iw)
            if (dummy.eq.-1.0) then
c              write(6,*) irec,nr,nmod,dummy,iw,' lambda vector scaling'
              do k = 1,iw
                 wave(k) = vec(k)
              enddo
            else
               n = n + 1
               age(n) = dummy        
               do k = 1, iw
                  if (vec(k).le.0) then
                       fl(n,k) = 0
                  else
                       fl(n,k) =  vec(k)
                  endif 
               enddo
            endif
         enddo
         jtot = n
c   extinction (iass=1 ==> 0 extinction)
         if (i.ge.modext(1) .and. i.le.modext(2)) then
           iassm = nebv   
         else
           iassm = 1
         endif 
         do ie = 1, iassm                 ! loop on extinction
c           zfmax=0
c           zfmin=10000
c           do k = 1,nzf
c              if (zf(k).ge.zfmax) zfmax=zf(k)
c              if (zf(k).le.zfmin) zfmin=zf(k)
c           enddo   
c          if (nrec(i).eq.1) then 
c             iform = 1
c          else
c             iform = nzf 
c          endif   
c          do k = 1, iform                 ! loop on Zform 
              if (nrec(i).eq.1) then
                 zform=20
              else
                 zform=zf(i)
              endif   
c  z step 
             if  (zmax.le.6) then
                 iz = idnint(zmax/dz) + 1
                 izl = iz
             else
                 izl = idnint(6./dz) + 1
                 iz  = idnint(6./dz) + 1 + idnint((zmax-6.)/dzup)
             endif
c
             do m = 1 , iz                ! loop on z
                if (m.le.izl) then
                   zoss = dz*(m-1)
                else
                   zoss = dz*(izl-1) + dzup*(m-izl)
                endif   
                if (zoss.ge.zform) then
                   eta=0
                   do if = 1,imag
                      mag(if)   = 99.
                      kcor(if)  = 0
                      ecor(if)  = 0
                      kecor(if) = 0
                      kapcor(if)= 0
                   enddo   
                   goto 35 
                endif   
                toss = timy(zoss,h0,om0,l0)-timy(zform,h0,om0,l0)
c take oldest age according to zform for flux at z=0
                if (m.eq.1) then         ! at z=0
                   z0=0.
                   tform = timy(z0,h0,om0,l0)-timy(zform,h0,om0,l0)
                   pass=0
                   if (nrec(i).eq.1) jtot=1
                   do ma = 1,jtot
                      if (age(ma).le.toss) then
                         agemax = ma
                      endif 
                   enddo
                   do n = 1, iw
                      flux0(n) =  fl(agemax,n)                  
                   enddo          
                endif
c take closest age according to zform for flux at z
                if (nrec(i).eq.1) then
c only one age -> one spectrum at any z 
                  do n = 1, iw
                    flux(n) =  fl(1,n)                  
                  enddo
                else
c loop on age to search for the lower and higher z   
                  do ma = 1,jtot-1          ! loop on Age
                     eta=age(ma)-toss
                     if ((eta*(age(ma+1)-toss)).le.0.) then
                       eta=tform-age(ma)
                       zs=zeta(eta,h0,om0,l0)    
                       eta=tform-age(ma+1)
                       zi=zeta(eta,h0,om0,l0)    
                       do n = 1, iw
                         fluxs(n) =  fl(ma,n)                  
                         fluxi(n) =  fl(ma+1,n)                  
                       enddo
                       eta=timy(zoss,h0,om0,l0)-timy(zform,h0,om0,l0)
                       eta=age(ma)
c                       write(*,'(I4,1x,10(E12.3,1x))')
c     >   i,zi,zoss,zs
                       goto 33
                     endif 
                  enddo
                endif  
 33             continue
c  apply extinction if any 
                val=ebv(ie)
                if (val.gt.0) then
                    iextinc=iext(1) 
                    do j = 1,iextinc
                      extinc(1,j)=exti(1,1,j)
                      extinc(2,j)=exti(1,2,j)
                    enddo
                   if (pass.eq.0) then
                  call scale_ext(wave,flux0,iw,extinc,iextinc,val,fext) 
                     do n = 1, iw
                       flux0(n)=fext(n)
                     enddo              
                   endif
                   if (nrec(i).eq.1) then
                  call scale_ext(wave,flux,iw,extinc,iextinc,val,fext) 
                     do n = 1, iw
                       flux(n)=fext(n)
                     enddo   
                   else
                  call scale_ext(wave,fluxs,iw,extinc,iextinc,val,fext) 
                     do n = 1, iw
                       fluxs(n)=fext(n)
                     enddo   
                  call scale_ext(wave,fluxi,iw,extinc,iextinc,val,fext) 
                     do n = 1, iw
                       fluxi(n)=fext(n)
                     enddo                      
                   endif                
                endif                
                call flush(6)
                write(6,1030) i,ebv(ie),eta,zoss,zform,char(13)
                call flush(6)
c
c computes mag(0,t(0)) -> mag1
                 if (pass.eq.0) then
                     z0=0.
                     call scale_opa(wave,flux0,iw,opal,opat,iopa,z0,
     >                     wnew,fnew,wnmax)
                     call  colgrid(fnew,wnew,wnmax,z0,
     >               flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag1)
                     pass=1
                 endif     
                 if (nrec(i).eq.1) then
c  computes mag(z,t(z)) -> mag
                  call scale_opa(wave,flux,iw,opal,opat,iopa,zoss,
     >                     wnew,fnew,wnmax)
                  call  colgrid(fnew,wnew,wnmax,zoss,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag)
c  computes mag(0,t(z)) -> mag2
                  z0=0.
                  call scale_opa(wave,flux,iw,opal,opat,iopa,z0,
     >                     wnew,fnew,wnmax)
                  call  colgrid(fnew,wnew,wnmax,z0,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag2)
c computes mag(z,t(0)) -> mag0
                  call scale_opa(wave,flux0,iw,opal,opat,iopa,zoss,
     >                     wnew,fnew,wnmax)
                  call  colgrid(fnew,wnew,wnmax,zoss,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag0)
c  Computes the various corrections 
                  z0=0.
                  distmod0= funz(z0,h0,om0,l0)
                  distmodz= funz(zoss,h0,om0,l0)
                  do if=1,imag
c  computes kcor(z)  = Mabs(z,t(0)) - Mabs(0,t(0))
                     kcor(if)   = mag0(if)-distmodz
     >                          -(mag1(if)-distmod0) 
c  computes ecor(z)  = Mabs(z,t(z)) - Mabs(z,t(0))
                     ecor(if)   = mag(if) - mag0(if)
c  computes kecor(z) = Mabs(z,t(0)) - Mabs(0,t(0))
                     kecor(if)  = kcor(if) + ecor(if) 
c  computes kapcor(z)= Mabs(z,t(z)) - Mabs(0,t(z))
                     kapcor(if) = mag(if)-distmodz
     >                          -(mag2(if)-distmod0) 
                  enddo                  
                 else                      ! nrec(i)>1
c mag(z,t(z)) -> magi,mags
                 call scale_opa(wave,fluxi,iw,opal,opat,iopa,zi,
     >                     wnew,fnew,wnmax)
                 call  colgrid(fnew,wnew,wnmax,zi,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,magi)
                 call scale_opa(wave,fluxs,iw,opal,opat,iopa,zs,
     >                     wnew,fnew,wnmax)
                 call  colgrid(fnew,wnew,wnmax,zs,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mags)
c mag(0,t(z)) -> mag2i,mag2s
                 z0=0.
                 call scale_opa(wave,fluxi,iw,opal,opat,iopa,z0,
     >                     wnew,fnew,wnmax)
                 call  colgrid(fnew,wnew,wnmax,z0,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag2i)
                 call scale_opa(wave,fluxs,iw,opal,opat,iopa,z0,
     >                     wnew,fnew,wnmax)
                 call  colgrid(fnew,wnew,wnmax,z0,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag2s)
c computes mag(z,t(0)) -> mag0i,mag0s
                 call scale_opa(wave,flux0,iw,opal,opat,iopa,zi,
     >                     wnew,fnew,wnmax)
                 call  colgrid(fnew,wnew,wnmax,zi,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag0i)
                 call scale_opa(wave,flux0,iw,opal,opat,iopa,zs,
     >                     wnew,fnew,wnmax)
                 call  colgrid(fnew,wnew,wnmax,zs,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag0s)
c  Computes the various corrections 
                 distmod0 = funz(z0,h0,om0,l0)
                 distmodzi= funz(zi,h0,om0,l0)
                 distmodzs= funz(zs,h0,om0,l0)
                 do if=1,imag
c  computes kcor(z)  = Mabs(z,t(0)) - Mabs(0,t(0))
                    kcori   = mag0i(if)-distmodzi
     >                           -(mag1(if)-distmod0)
                    kcors   = mag0s(if)-distmodzs
     >                           -(mag1(if)-distmod0)
                    kcor(if)= interp(zi,kcori,zs,kcors,zoss)
c  computes ecor(z)  = Mabs(z,t(z)) - Mabs(z,t(0))
                    ecori   = magi(if) - mag0i(if)
                    ecors   = mags(if) - mag0s(if)
                    ecor(if)=  interp(zi,ecori,zs,ecors,zoss)
c  computes kecor(z) = Mabs(z,t(0)) - Mabs(0,t(0))
                    kecor(if)  = kcor(if) + ecor(if) 
c  computes kapcor(z)= Mabs(z,t(z)) - Mabs(0,t(z))
                    kapcori = magi(if)-distmodzi
     >                           -(mag2i(if)-distmod0)
                    kapcors = mags(if)-distmodzs
     >                           -(mag2s(if)-distmod0)
                    kapcor(if) = interp(zi,kapcori,zs,kapcors,zoss)
c  computes mag 
                     kcori  = magi(if)
                     kcors  = mags(if)
                    mag(if) = interp(zi,kcori,zs,kcors,zoss)
                 enddo                     
                 endif           ! end if nrec(i)>1 
c                     
c   writing in output file 
 35              continue 
                 do if = 1,imag
                   mag(if) = mag(if) -2.5*dlog10(fcorr(if))
                 enddo
                 orec = orec + 1
         write(2,rec=orec) nmod,ebv(ie),zoss,eta,zform,imag,
     >                (mag(if),if=1,imag),(kecor(if),if=1,imag)
         if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y") then
            write(5,1010) nmod,ebv(ie),zoss,eta,zform,imag
     >  ,(mag(if),if=1,imag),(kcor(if),if=1,imag)
     >  ,(ecor(if),if=1,imag),(kecor(if),if=1,imag)
     >  ,(kapcor(if),if=1,imag)      
         endif   
             enddo               ! end loop on  z
c          enddo                  ! end loop on  zform
         enddo                   ! end loop on extinction 
      enddo                      ! end loop on SEDs 
c
      write(6,*) "            "
      write(6,*) " DONE       "
c    wrinting INFOs in doc file ...
      cr_date=fdate()
      write(4,'(2A)')     "CONFIG_FILE      ",config(1:lnblnk(config))
      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then   
        write(4,'(2A)')     "LIB_TYPE         ","GALAXY"
      elseif (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q") then
        write(4,'(2A)')     "LIB_TYPE         ","QSO"
      endif   
      write(4,'(2A)')      "LIB_NAME         ",lib(1:lnblnk(lib))
      write(4,'(A,3x,I8)') "NUMBER_SED       ",imas
      write(4,'(A,3x,I8)') "NUMBER_ROWS      ",orec
      write(4,'(A,3x,I8)') "RECORD_LENGTH    ",orecmax    
      write(4,'(101(A,2x))')"FILTER_FILE      ",
     >    filters(1:lnblnk(filters))
      write(4,'(101(A,2x))')"FILTERS          ",
     >                           (name(j)(1:lnblnk(name(j))),j=1,imag)
      write(4,'(51A)')     "MAG_TYPE         ",magtyp(1:5)
      write(4,'(A,100(f8.3,1x))')
     >                     "AB_COR           ",(abcor(j),j=1,imag)
      write(4,'(A,3(f7.3,1x))')"Z_STEP           ",dz,zmax,dzup   
      write(4,'(A,100(f6.3,1x))')"Z_FORM           ",(zf(if),if=1,nzf)
      write(4,'(A,3(f7.3,1x))')"COSMOLOGY        ",h0,om0,l0
      write(4,'(2A)')      "EXTINC_LAW       ",extlaw(1:lnblnk(extlaw))
      write(4,'(A,100(f6.3,1x))')"EB_V             ",(ebv(i),i=1,nebv)
      write(4,'(A,2(I6,1x))')  "MOD_EXTINC       ",modext(1),modext(2)
      write(4,'(2A)')     "CREATION_DATE    ",cr_date(1:lnblnk(cr_date))
      close(1)
      close(2)
      close(4)
      if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y") close(5)
c 999  format("# EB_V        : ",10(f8.3,1x))
 1010 format(i6,1x,f8.3,1x,f8.3,1x,E12.6,1x,f8.3,1x,i3,1x,500(f10.4,1x))
 1030 format("mod ->",I6," extinction:",f6.3," age(Gyr):",E12.6,
     >  " z: ",f7.3," zform: ",f7.3,a,$)
c
      stop
      end
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine colgrid(flux,wave,iw,zoss,
     >     flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag)
c
      implicit none
c
      integer*4 i,j,n,iw,imag,wmax,amax,maxsize,nbf
      INCLUDE 'dim_filt.decl'
c      parameter  (wmax=8000, amax=500,maxsize=110000)
      parameter  (amax=500,maxsize=110000)
      INCLUDE 'dim_wave.decl'      
      integer*4 nlrep,nlveg,jmax(nbf)
c
      real*8 mag(nbf),h0,om0,l0,zoss,z1,funz,fz
      real*8 rep(maxsize),lamb(maxsize)
      real*8 wave(wmax),flux(wmax),f(wmax),w(wmax)
      real*8 magab,magveg,abveg,abcor(nbf)
      real*8 lbveg(wmax),flveg(wmax)
      real*8 lambf(nbf,maxsize),repf(nbf,maxsize)
      character*512  magtyp

      external  funz
      common /univ/ h0,om0,l0
c
c  put flux and lambda according to z , exti opa...
      z1 = 1+zoss     
      fz = funz(zoss,h0,om0,l0)
      do n = 1 , iw
         f(n) =  flux(n)/z1                 ! erg.s-1.cm-2.A-1
         if (wave(n).le.500) f(n) = 0
         w(n) = wave(n)*z1 
      enddo
c
      do i =  1, imag
         nlrep=jmax(i)
         do j = 1,nlrep
	   lamb(j) = lambf(i,j)
           rep(j) = repf(i,j)
 	 enddo
c    check if lambda star pass fully through the filter
         if (lamb(nlrep).le.w(iw)) then
            call cal_mag(f,w,iw,rep,lamb,nlrep,
     >      flveg,lbveg,nlveg,magtyp,magveg,magab,abveg)
            abcor(i) = abveg
            if (magtyp(1:1) .eq. 'A')   mag(i)= magab + fz 
            if (magtyp(1:1) .eq. 'V')   mag(i)= magveg + fz
          else
            mag(i)=99
            abcor(i) = 99
          endif 
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 FUNCTION interp(x1,y1,x2,y2,x)
      implicit none 
c
      real*8  x1,y1,x2,y2,x,a,b
c
      a = (y1-y2)/(x1-x2)
      b = y1 -a*x1
      interp = a*x + b 
c
      return
      end




