c   last modif : 19/10/00
      program mag_eros
c     1 option : mag_gal  
c                        -c zphot.para    : config file
c     
c     Construction des differents modeles d'evolution pour
c     la mesure des z-photometriques fonction de Omo,Oml,Ho des
c     bandes photometriques et de l'extinction
c
c     Note: This program use all the Zformation for each model
c           The redshift does not fit with the zphot program !!
c            --> cannot be used for zphot 
      implicit none
      integer*4  nbf,maxsize
      integer*4  wmax,amax
c      parameter  (wmax=8000, amax=500)
      parameter  (amax=500)
      INCLUDE 'dim_wave.decl'      
      parameter  (maxsize=110000)
      INCLUDE 'dim_filt.decl'
      integer*4  jmax(nbf),pass    
      integer*4  imag,iass,imas,i,j,k,n,ie,if,iw,nebv,nzf
      integer*4  m,iassm,iform,test,agemax
      integer*4  recmax,nrec(10000),reclmax,irec
      integer*4  nmod,nr,jtot,orec,orecmax,nlveg,lnblnk
      integer*4  paravi(2),modext(2),iext,iopa(81),wnmax
c
      real*8     h0,om0,l0,val,dummy,funz
      real*8     timy,eta,zoss,distmod0,distmodz
      real*8     wave(wmax),exti(2,wmax)
      real*8     vec(wmax),zf(100)
      real*8     age(amax),fl(amax,wmax),flux(wmax),fext(wmax)
      real*8     paravr(100),ebv(100),mag(nbf),flux0(wmax)
      real*8     ecor(nbf),kcor(nbf),kecor(nbf),kapcor(nbf)
      real*8     mag0(nbf),mag1(nbf),mag2(nbf)
      real*8     lbveg(wmax),flveg(wmax)
      real*8     lambf(nbf,maxsize), repf(nbf,maxsize)
      real*8     opal(81,wmax),opat(81,wmax)
      real*8     wnew(wmax),fnew(wmax)
      real*8     zp(nbf),flmoy(nbf),flwidth(nbf),abcor(nbf)
      real*8     zfmax,zfmin,zform,tform,toss,zeta,z0
      real*8     fcorr(nbf)
c
      character*512  name(nbf),magtyp,param,vega
      character*512 paravc(100),valc,fileop(100),sedtyp,str,zpwork
      character*512  file,filters,lib,colib,zpdir,extlaw,config,outasc
      character*512 cr_date,fdate
      external   timy,funz,zeta
c
      common /univ/ h0,om0,l0
c      common /pas/ dz
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
      if (test.ne.1 .or. extlaw(1:4).eq.'NONE') then
          if (test.ne.1 ) call err_option(param,2)
          extlaw='NONE'
          nebv=1
          ebv(1) =0
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
c
         param='-MOD_EXTINC'
         call geti_option(param,config,2,paravi,test)
         if (test.eq.2)  modext(1) = paravi(1)
         if (test.eq.2)  modext(2) = paravi(2)
         if (test.ne.2)  call err_option(param,2)
         if (test.ne.2)  modext(1) = 0
         if (test.ne.2)  modext(2) = 0
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

c      param='-Z_STEP'
c      call getf_option(param,config,3,paravr,test)
c      if (test.eq.3)  dz = paravr(1)
c      if (test.eq.3)  zmax = paravr(2)
c      if (test.eq.3)  dzup = paravr(3)
c      if (test.ne.3)  call err_option(param,2)
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
       write(6,'(A,100(f5.2,1x))') "# Z_FORM     : ",(zf(i),i=1,nzf) 
       write(6,'(A,3(f7.3,1x))') "# COSMOLOGY   : ",h0,om0,l0
       write(6,'(2A)') "# EXTINC_LAW  : ",extlaw(1:lnblnk(extlaw))   
       write(6,'(A,100(f6.3,1x))') "# EB_V        :",(ebv(i),i=1,nebv)
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
      write(6,*) 'number of record :',imas,recmax,reclmax
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
      if (extlaw(1:4) .ne. 'NONE' .and. extlaw(1:4) .ne. 'none') then
c  read extinction law 
        file = zpdir(1:lnblnk(zpdir)) //'/ext/' 
     >   // extlaw(1:lnblnk(extlaw))
        open(1,file=file)
        i=0
        do while (.true.)
          i=i+1
          read(1,*,end=11) exti(1,i),exti(2,i)
        enddo
 11     iext=i-1
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
      write(6,650) imag,(name(j)(1:lnblnk(name(j))),j=1,imag)
 650  format('* ',I2,' filters selected : ',100(A10,1x))
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
      orecmax = 8*(5+imag*2)   ! if kcorrection added 
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
        write(5,'(A,1x,I6,1x,A,1x,20A)') "#",imag,magtyp(1:5),
     > (name(j)(1:lnblnk(name(j))),j=1,imag)
      endif
cccccccccccccccccccccccccccccccccccccccccccccc
c   Computing Magnitudes        cccccccccccccc   
cccccccccccccccccccccccccccccccccccccccccccccc
      irec = 0
      orec = 0
      do i = 1, imas  ! number of SEDs 
c  flux and lambda for each SEd and number of age (n)
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
c
c  computing mag for each SED at all z (with step of dz up to 6
c   and 0.2 up to zmax if zmax>6
c  and extinction (iass=1 ==> 0 extinction)
         if (i.ge.modext(1) .and. i.le.modext(2)) then
           iassm = nebv   
         else
           iassm = 1
         endif 
         do ie = 1, iassm  ! extinction
           zfmax=0
           zfmin=10000
           do k = 1,nzf
              if (zf(k).ge.zfmax) zfmax=zf(k)
              if (zf(k).le.zfmin) zfmin=zf(k)
           enddo   
          if (nrec(i).eq.1) then 
             iform = 1
          else
             iform = nzf 
          endif   
          do k = 1, iform
             if (nrec(i).eq.1) then
                zform=zfmax
             else
                zform=zf(k)
             endif   
             z0 = 0.
             tform = timy(z0,h0,om0,l0)-timy(zform,h0,om0,l0)
             pass=0
c  check the oldest
             if (nrec(i).eq.1) jtot=1
             agemax=1
             do m = 1,jtot
               eta=age(m)
               if (eta.le.tform) then
                 agemax = m
               endif 
             enddo
             do n = 1, iw
                flux0(n) =  fl(agemax,n)                  
             enddo             
c change jtot for ASCII file -> 1 age to 20 ==> 100 step in z
             if (nrec(i).eq.1) then
                jtot = 100
                do m = 1, jtot
                   age(m) = age(agemax)
                enddo   
             endif   
             do m = 1,jtot    ! Age
               eta = age(m)
               if (eta.le.tform) then
                 if (nrec(i).gt.1) then  
c if defined with age 
                   toss = tform - eta
                   zoss = zeta(toss,h0,om0,l0)                 
                   do n = 1, iw
                     flux(n) =  fl(m,n)                  
                   enddo
                 else 
c if defined with 1 age 
                   zoss = zform*(1.-dble(m)/dble(jtot))
                   do n = 1, iw
                     flux(n) =  fl(agemax,n)                  
                   enddo                  
                 endif 
                val=ebv(ie)
                if (val.gt.0) then
                  call scale_ext(wave,flux,iw,exti,iext,val,fext) 
                  do n = 1, iw
                    flux(n)=fext(n)
                  enddo   
                  if (pass.eq.0) then
                    call scale_ext(wave,flux0,iw,exti,iext,val,fext) 
                    do n = 1, iw
                      flux0(n)=fext(n)
                    enddo              
                  endif
                endif                
                call flush(6)
                write(6,1030) i,ebv(ie),eta,zoss,zform,char(13)
                call flush(6)
c computing various mag and correction kcor , ecor, kecor, kapcor
c
c  computes mag(z,t(z)) -> mag
                  call scale_opa(wave,flux,iw,opal,opat,iopa,zoss,
     >                     wnew,fnew,wnmax)
                  call  colgrid(fnew,wnew,wnmax,zoss,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag)
c
c  computes mag(0,t(z)) -> mag2
                  z0=0.
                  call scale_opa(wave,flux,iw,opal,opat,iopa,z0,
     >                     wnew,fnew,wnmax)
                  call  colgrid(fnew,wnew,wnmax,z0,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag2)
c
c computes mag(z,t(0)) -> mag0
                  call scale_opa(wave,flux0,iw,opal,opat,iopa,zoss,
     >                     wnew,fnew,wnmax)
                  call  colgrid(fnew,wnew,wnmax,zoss,
     >              flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag0)
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
c
c
c  Computes the various corrections 
                     z0=0.
                     distmod0= funz(z0,h0,om0,l0)
                     distmodz= funz(zoss,h0,om0,l0)
                   do if=1,imag
c  computes kcor(z)  = Mabs(z,t(0)) - Mabs(0,t(0))
                    kcor(if)   = mag0(if)-distmodz
     >                         -(mag1(if)-distmod0) 
c  computes ecor(z)  = Mabs(z,t(z)) - Mabs(z,t(0))
                    ecor(if)   = mag(if) - mag0(if)
c  computes kecor(z) = Mabs(z,t(0)) - Mabs(0,t(0))
                    kecor(if)  = kcor(if) + ecor(if) 
c  computes kapcor(z)= Mabs(z,t(z)) - Mabs(0,t(z))
                    kapcor(if) = mag(if)-distmodz
     >                         -(mag2(if)-distmod0) 
                   enddo                  
c                     
c   writing in output file 
                 do if=1,imag
                    mag(if) = mag(if) -2.5*dlog10(fcorr(if))
                 enddo
                 orec = orec + 1
         write(2,rec=orec) nmod,ebv(ie),zoss,eta,imag,
     >                (mag(if),if=1,imag),(kcor(if),if=1,imag)
         if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y") then
            write(5,1010) nmod,ebv(ie),zoss,eta,imag,zform,
     >                (mag(if),if=1,imag),(kcor(if),if=1,imag)
     >  ,(ecor(if),if=1,imag),(kecor(if),if=1,imag)
     >  ,(kapcor(if),if=1,imag)      
                 endif   
c                endif
               endif  
             enddo
           enddo
         enddo   
 35      continue
      enddo
c
      write(6,*) "            "
      write(6,*) " DONE       "
c   wrinting INFOs in doc file ...
      cr_date=fdate()
      write(4,'(2A)')     "CONFIG_FILE      ",config(1:lnblnk(config))
c
      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then   
        write(4,'(2A)')     "LIB_TYPE         ","GALAXY"
      elseif (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q") then
        write(4,'(2A)')     "LIB_TYPE         ","QSO"
      endif   
      write(4,'(2A)')      "LIB_NAME         ",lib(1:lnblnk(lib))
      write(4,'(A,3x,I8)') "NUMBER_SED       ",imas
      write(4,'(A,3x,I8)') "NUMBER_ROWS      ",orec
      write(4,'(A,3x,I8)') "RECORD_LENGTH    ",orecmax    
      write(4,'(51(A,2x))')"FILTER_FILE      ",
     >    filters(1:lnblnk(filters))
      write(4,'(51(A,2x))')"FILTERS          ",
     >                           (name(j)(1:lnblnk(name(j))),j=1,imag)
      write(4,'(51A)')     "MAG_TYPE         ",magtyp(1:5)
      write(4,'(A,100(f8.3,1x))')
     >                     "AB_COR           ",(abcor(j),j=1,imag)
      write(4,'(A,100(f7.2,1x))')"Z_FORMATION      ",(zf(if),if=1,nzf)    
      write(4,'(A,3(f7.3,1x))')"COSMOLOGY        ",h0,om0,l0
      write(4,'(2A)')     "EXTINC_LAW       ",extlaw(1:lnblnk(extlaw))   
      write(4,'(A,100(f6.3,1x))')"EB_V             ",(ebv(i),i=1,nebv)
      write(4,'(A,2(I6,1x))')  "MOD_EXTINC       ",modext(1),modext(2)
      write(4,'(2A)')     "CREATION_DATE    ",cr_date(1:lnblnk(cr_date))
c
      close(1)
      close(2)
      close(4)

      if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y") close(5)
 999  format("# EB_V        : ",10(f8.3,1x))
 1010 format(i6,1x,f5.3,1x,f6.3,1x,E12.6,1x,i3,1x,500(f10.4,1x))
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
      INCLUDE 'dim_wave.decl'      
c      parameter  (wmax=8000, amax=500,maxsize=110000)
      parameter  (amax=500,maxsize=110000)
      integer*4 nlrep,nlveg,jmax(nbf)
c
      real*8 mag(nbf),h0,om0,l0,zoss,z1,funz,fz
      real*8 rep(maxsize),lamb(maxsize)
      real*8 wave(wmax),flux(wmax),f(wmax),w(wmax)
      real*8 magab,magveg,abveg,abcor(nbf)
      real*8 lbveg(wmax),flveg(wmax)
      real*8 lambf(nbf,maxsize),repf(nbf,maxsize)
      character*512     magtyp

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




