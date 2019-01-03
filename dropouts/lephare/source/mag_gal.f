c   last modif : 19/10/00
      program mag_gal
c     2 options : mag_gal  
c                        -c zphot.para    : config file
c                        -t (G/g or Q/q)  : GAL/QSO library
c     
c     Construction des differents modeles d'evolution pour
c     la mesure des z-photometriques fonction de Omo,Oml,Ho des
c     bandes photometriques et de l'extinction
c
      implicit none
      integer*4  nbf,maxsize
      integer*4  wmax,amax
c      parameter  (wmax=8000, amax=500)
      parameter  (amax=500)
      INCLUDE 'dim_wave.decl'      
      parameter  (maxsize=110000)
      INCLUDE 'dim_filt.decl'
      INCLUDE 'out_unit.decl' 
      integer*4  jmax(nbf)    
      integer*4  imag,iass,imas,i,j,k,n,ifilt,iw
      integer*4  iform,m,iforml,test
      integer*4  recmax,nrec(maxsize),reclmax,irec
      integer*4  nmod,nr,jtot,orec,orecmax,nlveg,lnblnk
      integer*4  paravi(500),iopa(81),wnmax,librec(amax)
c
      real*8     h0,om0,l0,dz,zmax,dzup,dummy,funz,dist_mod
      real*8     timy,tuniv,eta,zoss,valebv
      real*8     wave(wmax)
      real*8     vec(wmax)    ! ,z0
      real*8     age(amax),fl(amax,wmax),flux(wmax),fext(wmax),tau(amax)
      real*8     wir(wmax),hir0(wmax),hir(wmax),ltir0,ltir1,ltir(amax)
      real*8     paravr(500),mag(nbf)
      real*8     magko(nbf),kcor(nbf)
      real*8     lbveg(wmax),flveg(wmax)
      real*8     lambf(nbf,wmax), repf(nbf,wmax)
      real*8     opal(81,wmax),opat(81,wmax)
      real*8     wnew(wmax),fnew(wmax)
      real*8     zp(nbf),flmoy(nbf),flwidth(nbf),abcor(nbf)
      real*8     fcorr(nbf)
c
      integer*4   pass
      real*8      zrec(3,1000)
c
      character*4096 name(nbf),magtyp,param
      character*4096 paravc(500)
      character*4096 fileop(500),sedtyp,zpwork
      character*4096 file,filters,lib,colib,zpdir,config,outasc
      character*4096 cr_date,fdate
      character*4096 str,libphys
c extinction
      character*4096 extlaw(10)
      integer*4      modext(20),nmodext,nextlaw,nebv,iext(10)
      integer*4      iextinc,iel,iebv
      real*8         ebv(500),exti(10,2,wmax),extinc(2,wmax)
c physical parameters 
      integer*4      kphys,recp,modp
      real*8         physpara(50,maxsize),nuvr(amax),NUV_R
      real*8         magkos(maxsize,nbf),ldust(amax)
c Emission lines
      real*8         em(nbf),eml(nbf),NUVabs
      real*8         ext_em(7)
      character*4096 emlines
      integer*4      imuv 
      real*8         lambf_UV(wmax),repf_UV(wmax),abcor_UV

      external      timy,funz
      real*8      Lsol,pc,c
      parameter  (Lsol=3.826e33,pc=3.086e18)
      parameter  (c=2.99792458e+18)      
     
      common /univ/ h0,om0,l0
      common /pas/ dz
      common /input/ iass, imas

c      If you want to give a name to the output screen file
      if(UO.eq.30)
     .      open(UO,file='/tmp/screenMagGal.dat',status='unknown') 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc   PARAMETERS                         cccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  environmental variable 
      call getenv('LEPHAREDIR',zpdir)
      test=lnblnk(zpdir)
      if (test .eq. 0) then
        write(UO,*) 'WARNING :  variable LEPHAREDIR not defined'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(UO,*) 'WARNING :  variable LEPHAREWORK not defined'
        stop
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c help on line
      param='mag_gal'
      call get_help(test)
      if (test .eq. 1) call help(param)
cccccccccccccccccccccccccccccccccccccccccccccccc
c Initialisation of Input parameter from config file
      param='-c'
      call get_conf(param,config,test)
      if (test.ne.1) call err_option(param,1)
      call get_path(config)
      param='-t'
      call get_conf(param,sedtyp,test)
      if (test.ne.1) call err_option(param,1)
      if (sedtyp(1:1) .eq.'g') sedtyp = "GAL"
      if (sedtyp(1:1) .eq.'q') sedtyp = "QSO"
      if (sedtyp(1:1) .eq.'s') sedtyp = "STAR"
c  read option
      param='-COSMOLOGY'
      call getf_option(param,config,3,paravr,test)
      if (test.eq.3)   h0=paravr(1)
      if (test.eq.3)   om0=paravr(2)
      if (test.eq.3)   l0=paravr(3)
      if (test.ne.3)  call err_option(param,1)

      param='-FILTER_FILE'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  filters=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)

      param='-MAGTYPE'
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1)  magtyp=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)

c   Max of extinction law  is 10
      nextlaw=0 
      do i = 1, 10 
        extlaw(i) = 'NONE'
      enddo
      nebv=1
      do i = 1, 500
        ebv(i) = 0. 
      enddo
      nmodext=0
      do i = 1, 20
         modext(i) = 0 
      enddo
c
      param='-EXTINC_LAW'
      call  getc_option(param,config,500,paravc,test)
c      write(6,*)  'extinction ',test,'  ',
c     >             (paravc(i)(1:lnblnk(paravc(i))+1),i=1,10)
      if (test.eq.0) then
c         call err_option(param,2)
         extlaw(1)='NONE'
         nextlaw=0
      else
        if (test.gt.10) test=10
        nextlaw=test
        do  i = 1, test 
            extlaw(i)=paravc(i)(1:lnblnk(paravc(i)))
        enddo
      endif 
      if ( extlaw(1)(1:lnblnk(extlaw(1))) .eq. 'NONE') nextlaw=0 
c  Read E(B-V) and models with_extinction
      if (nextlaw .eq. 0 ) then 
         nebv=1
         ebv(1)=0.
         modext(1) = 0
         modext(2) = 0 
         nmodext =0
      else   
         param='-EB_V'
         call getf_option(param,config,500,paravr,test)
         if (test.ge.1) then
           nebv=test
           do i = 1,nebv 
              ebv(i) =  paravr(i)
           enddo
         else
            nebv=1
            ebv(1)=0.
            extlaw(1)='NONE'
            nextlaw=0
            modext(1) = 0
            modext(2) = 0 
            nmodext =0            
         endif   
c  read model ranges with extincion 
         param='-MOD_EXTINC'
         call geti_option(param,config,500,paravi,test)         
         if (test.eq.0 .OR. test .ne. (2*nextlaw)  ) then 
            modext(1) = 0
            modext(2) = 0
            nmodext=0
            nebv=1
            ebv(1)=0.
            extlaw(1)='NONE'
            nextlaw=0
c           
            write(UO,*) ' '
            write(UO,*) ' No extinction applied '
            if (test .ne. (2*nextlaw) ) 
     > write(UO,*) ' MOD_EXTINC  not consistent with number of  
     > extinction law in EXTINC_LAW '   
            write(UO,*) ' '
         else
            nmodext=test
            do i = 1, test 
               modext(i) = paravi(i)
            enddo
         endif
      endif
c  Z
      param='-Z_STEP'
      call getf_option(param,config,3,paravr,test)
      if (test.eq.3)  dz = paravr(1)
      if (test.eq.3)  zmax = paravr(2)
      if (test.eq.3)  dzup = paravr(3)
      if (test.ne.3)  call err_option(param,1)
c library 
      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
        param='-GAL_LIB_IN'
      elseif (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q") then
        param='-QSO_LIB_IN'
      endif 
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1)  lib=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
         param='-GAL_LIB_OUT'
      elseif (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q") then
         param='-QSO_LIB_OUT'
      endif
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1)  colib=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
      param='-LIB_ASCII'
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1)  outasc=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,2)
      if (test.ne.1)  outasc='NO'
c emission lines
      emlines="NO" 
      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
         param='-EM_LINES'
         call  getc_option(param,config,1,paravc,test)
         if (test .eq. 1) then
            if ( paravc(1)(1:1) .eq. "y" .or.
     >           paravc(1)(1:1) .eq. "Y")  emlines= "YES"
         endif
         if (test .ne. 1) call err_option(param,2)
      endif        
cccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccc
c      READING FILES  
ccccccccccccccccccccccccccccccc
c  read doc file from INPUT LIBRARY 
      write(UO,*) ' reading library doc ...'
c
      file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     > // lib(1:lnblnk(lib)) // '.doc'
c
      param='LIB_FMT'
      call read_para(param,file,paravc,test)
      libphys='NO'
      if (test.eq.1) then 
c   setup automatically some  options for BC_STOCH 
         if(paravc(1)(1:lnblnk(paravc(1))).eq.'BC_STOCH') then
            libphys= 'YES'
            emlines= 'NO'
            outasc = 'NO'
            nebv=1
            ebv(1)=0.
            extlaw(1)='NONE'
            nextlaw=0
            modext(1) = 0
            modext(2) = 0 
            nmodext =0            
            write(UO,*) ' Automatic setup for LIB_FMT=BC_STOCH library'
            write(UO,*) '-EM_LINES=NO -EXTINC_LAW=NONE -LIB_ASCII=NO'
         endif
c   setup automatically some  options for FIR library 
         if(paravc(1)(1:lnblnk(paravc(1))).eq.'LW') then
            emlines= 'NO'
            nebv=1
            ebv(1)=0.
            extlaw(1)='NONE'
            nextlaw=0
            modext(1) = 0
            modext(2) = 0 
            nmodext =0            
            write(UO,*) ' Automatic setup for LIB_FMT=LW (FIR library)'
            write(UO,*) '-EM_LINES=NO  -EXTINC_LAW=NONE'
         endif

      else
         write(UO,*) ' LIB_FMT keyword not available in',
     >               file(1:lnblnk(file)),' Old library -> STOP'
         STOP
      endif
      param='NUMBER_SED'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') imas
c      write(UO,*) imas,'  ',libphys(1:2)
      if (libphys(1:2) .eq. 'NO') then 
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
      else
         do i = 1, imas
             nrec(i)=1
         enddo    
      endif
      param='NUMBER_ROWS'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') recmax
      param='RECORD_LENGTH'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') reclmax
c
      write(UO,*) 'number of record :',imas,recmax,reclmax
c  read Vega spectrum
      file = zpdir(1:lnblnk(zpdir)) // '/vega/VegaLCB.sed'
      open(1, file=file,status='old',err=56)
      read(1,*)   
      i=0
      do while (.true.)
	i=i+1
        read(1,*,end=10)  lbveg(i), flveg(i)
      enddo
 10   nlveg=i-1  
      close(1)
c  read NUV filter for Emission line 
      if (emlines(1:1) .eq. "Y") then 
         file = zpdir(1:lnblnk(zpdir)) //'/filt/galex/NUV.pb'
         open(1,file=file,status='old',err=56)
         read(1,*) 
         j=0
         do while (.true.)
          j=j+1
          read(1,*,end=14)  lambf_UV(j),repf_UV(j)
         enddo
 14      imuv=j-1  
         close(1)
         abcor_UV=1.76
      endif  
c read physical parameters from gallib in lib_bin
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if((sedtyp(1:1).eq."G" .or. sedtyp(1:1).eq."g") ) then 
           file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >           // lib(1:lnblnk(lib)) // '.phys'
          write(UO,*) 'reading Physical parameters from',
     >    file(1:lnblnk(file)),' with',recmax         
          if (recmax .gt. 0 ) then 
              open(1,file=file,status='unknown',err=56)
              k = 0
              do while (.true.)
                 k = k + 1
                read(1,*,end=15) recp,modp,(physpara(i,k),i=1,10)
              enddo  
 15           kphys = k - 1
             close(1)
             if (kphys.ne.recmax) 
     >       write(UO,*) ' Size problem for the  lib_bin library ',
     >                    lib(1:lnblnk(lib)),kphys,recmax          
          endif
      endif 
c
c  read the filters from  file
      file = zpwork(1:lnblnk(zpwork)) //'/filt/'
     >        // filters(1:lnblnk(filters))
      open(1,file=file,status='old',err=56)
      read(1,'(A)') str 
      call val_string(str,paravc,test)
      read(paravc(2),'(i10)') imag
      do i =  1, imag
         read(1,'(A)') str
         call val_string(str,paravc,test)
         read(paravc(2),'(i10)') jmax(i)
         name(i)=paravc(3)
         do j = 1,jmax(i)
	    read(1,*)  lambf(i,j), repf(i,j)
 	 enddo
      enddo
      close(1) 
c  getting AB corrections for info 
      call zeropoint(filters,zp,abcor,flmoy,flwidth,fcorr,imag)
c  read extinction laws 
      if (nextlaw .gt. 0 .and.  nmodext .gt. 0 ) then 
        do j = 1,nextlaw
           file = zpdir(1:lnblnk(zpdir)) //'/ext/' 
     >     // extlaw(j)(1:lnblnk(extlaw(j)))
           open(1,file=file,status='old',err=56)
           i=0
           do while (.true.)
             i=i+1
             read(1,*,end=11) exti(j,1,i),exti(j,2,i)
           enddo
 11        iext(j)=i-1
           close(1)
c          write(6,*) 'iext(j)=',iext(j) 
        enddo
      endif  

c  reading file with extragalacitic opacity
      file = zpdir(1:lnblnk(zpdir)) // '/opa/OPACITY.dat'
      open(1,file=file,status='old',err=56)
      do i = 1,81
         read(1,'(A)') str
         call val_string(str,paravc,test)
         fileop(i)=paravc(2)
      enddo
      close(1)
      do i = 1, 81
          file = zpdir(1:lnblnk(zpdir)) // '/opa/' 
     > // fileop(i)
          open(1,file=file,status='old',err=56)
          k = 0
          do while (.true.)
             k = k + 1
             read(1,*,end=12) opal(i,k), opat(i,k)
          enddo  
 12       iopa(i) = k - 1
          close(1)
      enddo    
      write(UO,650) imag,(name(j),j=1,imag)
 650  format('* ',I2,' filters selected : ',100(A10,1x))
      write(UO,900) imas
 900  format('* number of  models used =',i6)
c
cccccccccccccccccccccccccccccccccccccccccccccc
c  writing options    
       write(UO,'(A)')"############################################"
       write(UO,'(A)')"# It s computing the SYNTHETIC MAGNITUDES  #"
       write(UO,'(A)')"# For Gal/QSO libraries with these OPTIONS #"
       file = zpwork(1:lnblnk(zpwork)) //'/filt/'
     >      // filters(1:lnblnk(filters))
       write(UO,'(2A)')"# FILTER_FILE : ",file(1:lnblnk(file))  
       write(UO,'(2A)')"# MAGTYPE     : ",magtyp(1:lnblnk(magtyp))
       write(UO,'(A,50(f8.3,1x))')
     >                 "# FLUX_COR    : ",(fcorr(j),j=1,imag)
       if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
           file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >      // lib(1:lnblnk(lib)) //' (.doc & .bin & .phys)'
         write(UO,'(2A)') "# GAL_LIB_IN  : ",file(1:lnblnk(file)) 
         file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >      // colib(1:lnblnk(colib)) //'(.doc & .bin)'
         write(UO,'(2A)') "# GAL_LIB_OUT : ",file(1:lnblnk(file))     
       elseif (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q") then
           file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >      // lib(1:lnblnk(lib)) //'(.doc & .bin)'
         write(UO,'(2A)') "# QSO_LIB_IN  : ",file(1:lnblnk(file)) 
          file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >      // colib(1:lnblnk(colib)) //'(.doc & .bin)'
         write(UO,'(2A)') "# QSO_LIB_OUT : ",file(1:lnblnk(file))     
       endif  
       write(UO,'(A,3(f7.3,1x))') "# Z_STEP      : ",dz,zmax,dzup    
       write(UO,'(A,3(f7.3,1x))') "# COSMOLOGY   : ",h0,om0,l0
         
       if (nextlaw.gt.0)  write(UO,'(11(A,2x))') "# EXTINC_LAW  : ",
     >                    (extlaw(i)(1:lnblnk(extlaw(i))),i=1,nextlaw)
       if (nmodext .gt. 0) then 
            write(UO,'(A,20(I6,1x))') "# MOD_EXTINC  : ",
     >                                    (modext(i),i=1,nmodext)
       else 
           write(UO,'(A,20(I6,1x))') "# MOD_EXTINC  : ",
     >                                    (modext(i),i=1,2)
       endif 
       write(UO,'(A,10(f6.3,1x))') "# EB_V        :",(ebv(i),i=1,nebv)
       write(UO,'(2A)') "# EM_LINES    : ",emlines(1:lnblnk(emlines))
       write(UO,'(2A)') "# LIB_ASCII   : ",outasc(1:lnblnk(outasc)) 
       write(UO,'(A)')"############################################"
cccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c    OPENING INPUT and OUTPUT FILES           ccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Open input library with models (1) ...
      file = zpwork(1:lnblnk(zpwork)) // '/lib_bin/'
     >       // lib(1:lnblnk(lib)) // '.bin'
      open(1,file=file,status='unknown',access='direct',
     >                               recl=reclmax,err=56)
c  Open output mag lib file (2)
      if (emlines(1:1) .eq. "Y") orecmax=8*(8+imag*3)
      if (emlines(1:1) .eq. "N") orecmax=8*(8+imag*2)
      if (libphys(1:1) .eq. 'Y') orecmax=8*(7+imag*2)

c      orecmax = 4*(5+imag)
      file = zpwork(1:lnblnk(zpwork)) // '/lib_mag/'
     >    // colib(1:lnblnk(colib)) // '.bin'
      open(2,file=file,form='unformatted',status='unknown',
     > access='direct',recl=orecmax)
c  Open output doc file (4)
      file = zpwork(1:lnblnk(zpwork)) // '/lib_mag/'
     >    //colib(1:lnblnk(colib)) // '.doc'
      open(4,file=file,status='unknown')
c Open output mag ASCII file (5)
      if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y") then 
        file = colib(1:lnblnk(colib)) // '.dat'
        open(5,file=file,status='unknown')
        write(5,'(A,1x,I6,1x,A,1x,100A)') "#",imag,magtyp(1:5),
     > (name(j)(1:lnblnk(name(j)))," ",j=1,imag)
      endif
c
      NUV_R=10
      iform =1
      if (libphys(1:2) .eq. 'NO') then
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
               librec(n) = nr     
               do k = 1, iw
                  if (vec(k).le.0) then
                       fl(n,k) = 0
                  else
                       fl(n,k) =  vec(k)
                  endif 
               enddo  
               if (sedtyp(1:1).eq."G") then 
                  nuvr(n)  =-2.5*(physpara(2,nr)-physpara(3,nr))  
                  ldust(n) = physpara(5,nr)
                  ltir(n)  = ldust(n)
                  tau(n)   = physpara(9,nr)
c                  write(UO,*) nmod,' NUV-R= ',nuvr(n),age(n),ldust(n)
               endif
             endif
         enddo
c  number of ages per models  
         jtot = n
c  redshift info
         if (zmax.le.6.0) then
            iform = idnint(zmax/dz) + 1
            iforml = iform
         else
            iforml = idnint(6./dz) + 1
            iform = idnint(6./dz)  + 1 + idnint((zmax-6.)/dzup)
         endif
c  extinction 
         do iel = 0, nextlaw  ! extinction law
            if (iel.gt.0 .and. nmodext.eq.0) goto  35
            if (iel.gt.0 .and. nmodext.gt.0 .and.
     >     (i.lt.modext(2*iel-1) .or. i.gt.modext(2*iel)) ) goto 34 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c  check values for attenuation in Emis. lines
            if (emlines(1:1) .eq. "Y") then 
              if ( iel.ge.1) then 
                 iextinc=iext(iel)
                 do j = 1,iextinc
                    extinc(1,j)=exti(iel,1,j)
                    extinc(2,j)=exti(iel,2,j)
                 enddo
                 call ext_emlines(extinc,iextinc,
     >              lambf_UV,repf_UV,imuv,ext_em)
              else
                 do ifilt=1,7
                    ext_em(ifilt)= 0
                 enddo
              endif
            endif             
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
            do iebv = 1, nebv  ! ebv
              valebv=ebv(iebv)
              if (iel.eq.0 .and. iebv.gt.1 )    goto  34  ! do not apply ebv if extlaw=NONE
              if (iel.gt.0 .and. iebv.eq.1 )    goto  33  ! do not apply ebv=0 for extlaw
              do m = 1,jtot    ! age
                 if (ldust(m).gt.-10) then 
                   ltir(m)=ldust(m)
                 else
                   ltir(m) =-99 
                 endif
                 if (tau(m).gt.0 .and. age(m).gt.0) then 
                    if (age(m)/tau(m).gt.4 .and. valebv.ge.0.15) goto 33
                 endif
c  extract Restframe NUV-r with no extinction !
                 if (sedtyp(1:1).eq."G")  NUV_R=nuvr(m) 
c  extract flux at a given age 
                 do n = 1, iw
                    flux(n) =  fl(m,n)   
                 enddo
                 if (valebv .gt. 0) then
                    iextinc=iext(iel) 
                    do j = 1,iextinc
                      extinc(1,j)=exti(iel,1,j)
                      extinc(2,j)=exti(iel,2,j)
                    enddo
                    call scale_ext(wave,flux,iw,extinc,iextinc,
     >                                              valebv,fext) 
                    do n = 1, iw
                       flux(n)=fext(n)
                    enddo
                                      
c            compute Lir according to the different ext law if ltir does not exist in Phys para 
                    if (ltir(m).le.-90) then 
                      do n = 1, iw 
                        wir(n) = wave(n)/10000.        ! lbda in micron 
                        hir0(n) = fl(m,n)/(3.197e-7)   ! spectrum dust free
                        hir(n)  = flux(n)/(3.197e-7)   ! spectrum with extinction. 
c                                                      ! converting {erg/s/A/cm^2} to Lo/A  
                      enddo 
                      call TRAPZD(wir,hir,iw,wir(1),wir(iw),ltir1)
                      call TRAPZD(wir,hir0,iw,wir(1),wir(iw),ltir0)
               if (ltir0.gt.ltir1) ltir(m)=dlog10((ltir0-ltir1)*10000.)
c
                   endif
c                    
c   compute Attenuation in Em lines if Emission lines =Y 
c    dependance with extlaw and ebv not galaxy age : done only for first age
                 endif      
                 eta = age(m)*1.e-9         
                 call flush(6)
                 write(UO,1030) i,iel,valebv,eta,char(13)
                 call flush(6)
                 eta=age(m)
                 do k = 1, iform   ! z 
                    if (k.le.iforml) then
                       zoss = dz*(k-1)
                    else
                       zoss = dz*(iforml-1) + dzup*(k-iforml)
                    endif
                    tuniv = timy(zoss,h0,om0,l0)
                    dist_mod=funz(zoss,h0,om0,l0)
c                    if (eta.le.tuniv) then
                      call scale_opa(wave,flux,iw,opal,opat,iopa,zoss,
     >                     wnew,fnew,wnmax)
                      call  colgrid(fnew,wnew,wnmax,zoss,dist_mod,
     >                     flveg,lbveg,nlveg,lambf,repf,jmax,
     >                     imag,magtyp,mag)
c   computing k-correction 
                      if (zoss.le.1.e-5) then
                          do ifilt=1,imag
                             magko(ifilt) = mag(ifilt) 
                             kcor(ifilt)=0.
                          enddo
                      else
                          do ifilt=1,imag
                             kcor(ifilt)=mag(ifilt)-dist_mod
     >                                   -magko(ifilt)
                          enddo
                      endif
c   apply The correction factor for long wavelengths and various calib
                      do ifilt=1,imag
                      mag(ifilt) = mag(ifilt) -2.5*dlog10(fcorr(ifilt))
                      enddo
c   measure the flux from emission lines 
                      if (emlines(1:1) .eq. "Y" ) then 
                         do ifilt=1,imag
                           em(ifilt) =-999.d0
                           eml(ifilt) =-999.d0
                         enddo
c   only if object intrinsically blue (NUV-r)<4 (before extinction!)
                         if (NUV_R.le.4) then 
c                        write(*,*) nmod,iel,valebv,zoss,NUV_R
                            NUVabs=-20.d0 
c                            call addEmission(imag,NUVabs,zoss,valebv,
c     >                 dist_mod,flmoy,flwidth,lambf,repf,jmax,ext_em,em)
                            call addEmission2(imag,NUVabs,zoss,valebv,
     >                            dist_mod,lambf,repf,jmax,ext_em,em)
                            do ifilt=1,imag
                              if(em(ifilt).gt.0) then
                                eml(ifilt)=dlog10(em(ifilt))
                              else
                                eml(ifilt)=-999.
                              endif
                            enddo
                         endif 
                      endif
c  
c   writing output 
                      orec = orec + 1
                      if (emlines(1:1) .eq. "Y") then 
                 write(2,rec=orec) nmod,iel,valebv,ltir(m),zoss,
     >             eta,librec(m),imag,(mag(ifilt),ifilt=1,imag),
     >             (kcor(ifilt),ifilt=1,imag),(em(ifilt),ifilt=1,imag)
                      else
                 write(2,rec=orec) nmod,iel,valebv,ltir(m),zoss,
     >             eta,librec(m),imag,(mag(ifilt),ifilt=1,imag),
     >             (kcor(ifilt),ifilt=1,imag)
                      endif
c
                   if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y")then
                        if (emlines(1:1) .eq. "Y")
     >                    write(5,1020) nmod,iel,valebv,ltir(m),zoss,
     >            dist_mod,eta,librec(m),imag,(mag(ifilt),ifilt=1,imag),
     >           (kcor(ifilt),ifilt=1,imag),(eml(ifilt),ifilt=1,imag)
                        if (emlines(1:1) .eq. "N")
     >                    write(5,1010) nmod,iel,valebv,ltir(m),zoss,
     >            dist_mod,eta,librec(m),imag,(mag(ifilt),ifilt=1,imag),
     >            (kcor(ifilt),ifilt=1,imag)
                     endif
c
c                    endif  ! if age lower than Tuniv(z)
                  enddo    ! loop on z 
               enddo       ! loop on age 
 33            continue    ! do not repeat ebv=0 for different extinction law  
            enddo          ! loop on ebv  
 34         continue       ! do not apply this extinction law for this model 
         enddo             ! loop on extinction law    
 35      continue          ! do not apply extinction no model-range 
       enddo               ! loop on model 
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else                 ! For BC_STOCH   library 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       orec = 0
       read(1,rec=1) nr,nmod,dummy,iw,(vec(k),k=1,iw)
       if (dummy.eq.-1.0) then
         do k = 1,iw
            wave(k) = vec(k)
         enddo
         write(UO,*) nr,nmod,dummy,iw,' lambda vector extracted'
       endif
c
c       write(*,*) 'ok'
       if (zmax.le.6.0) then
          iform = idnint(zmax/dz) + 1
          iforml = iform
       else
          iforml = idnint(6./dz) + 1
          iform = idnint(6./dz)  + 1 + idnint((zmax-6.)/dzup)
       endif
c       write(*,*) iform,zmax,recmax
       do k = 1, iform                  ! z 
         if (k.le.iforml) then
            zoss = dz*(k-1)
         else
            zoss = dz*(iforml-1) + dzup*(k-iforml)
         endif
         zrec(1,k)=zoss
         if (k.gt.1) zrec(3,k-1)=orec
c         tuniv = timy(zoss,h0,om0,l0)
         dist_mod=funz(zoss,h0,om0,l0)
         pass=0
c         write(*,*) 'ok',recmax
         do i = 2, recmax              ! number of SEDs 
c            write(*,*) i,recmax
            read(1,rec=i) nr,nmod,dummy,iw,(vec(n),n=1,iw)
c            write(*,*) i,nr,nmod,dummy,iw,vec(1),vec(iw)
           if (dummy.eq.-1) then 
                write(UO,*) 'Pb with library -> STOP'
                STOP
            endif    
            eta = dummy
            iel=0
            valebv=0.d0 
            do n = 1, iw
               if (vec(n).le.0) then
                  flux(n) = 0
               else
                  flux(n) =  vec(n)
               endif 
            enddo  
c            write(UO,*) i,nr,nmod,dummy,iw,eta,flux(1),flux(iw)
c
c            if (eta.le.tuniv) then
               if (MOD(i,1000).eq.0) then  
                 call flush(6)
c               write(UO,'("Z= ",f6.3," model=",i8," age=",f9.6,a,$)')
c     >                 zoss,i,eta*1.e-9,char(13)
                 write(UO,'("Z= ",f6.3," model=",i8,a,$)')
     >                 zoss,i,char(13)
                 call flush(6)
               endif
c
               call scale_opa(wave,flux,iw,opal,opat,iopa,zoss,
     >                         wnew,fnew,wnmax)
c               write(*,*) 'ok2'
               call  colgrid(fnew,wnew,wnmax,zoss,dist_mod,
     >                     flveg,lbveg,nlveg,lambf,repf,jmax,
     >                     imag,magtyp,mag)
c               write(*,*) 'ok3'
c               do ifilt=1,imag
c                 mag(ifilt) = 20
c               enddo
               if (zoss.le.1.d-5) then
                  do ifilt=1,imag
                     magkos(i,ifilt) = mag(ifilt) 
                     kcor(ifilt)=0.
                  enddo
               else
                  do ifilt=1,imag
                     kcor(ifilt)=mag(ifilt)-dist_mod-magkos(i,ifilt)
                  enddo
               endif
               do ifilt=1,imag
                  mag(ifilt) = mag(ifilt) -2.5*dlog10(fcorr(ifilt))
               enddo
               orec = orec + 1
               write(2,rec=orec) nmod,iel,valebv,zoss,
     >             eta,nr,imag,(mag(ifilt),ifilt=1,imag),
     >            (kcor(ifilt),ifilt=1,imag)
c               write(*,*) i,orec,nmod,nr,zoss,eta,imag
               if (pass.eq.0) then
                     zrec(2,k)=orec
                     pass=1
               endif
c            endif   ! if age lower than Tuniv(z)
         enddo      ! loop on models 
         if (k.eq.iform) zrec(3,iform)=orec
        
       enddo         ! loop on z 
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(UO,*) "            "
      write(UO,*) " DONE       "
c   wrinting INFOs in doc file ...
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
      write(4,'(51(A,2x))')"FILTER_FILE      ",
     >    filters(1:lnblnk(filters))
      write(4,'(51(A,2x))')"FILTERS          ",
     >                           (name(j)(1:lnblnk(name(j))),j=1,imag)
      write(4,'(51A)')     "MAG_TYPE         ",magtyp(1:5)
      write(4,'(A,50(f8.3,1x))')
     >                     "AB_COR           ",(abcor(j),j=1,imag)
      write(4,'(A,50(f8.3,1x))')
     >                     "FLUX_COR         ",(fcorr(j),j=1,imag)
      write(4,'(A,3(f7.3,1x))')"Z_STEP           ",dz,zmax,dzup    
      write(4,'(A,3(f7.3,1x))')"COSMOLOGY        ",h0,om0,l0
      if (nextlaw.gt.0) then 
           write(4,'(11(A,2x))')      "EXTINC_LAW       ",
     >                    (extlaw(j)(1:lnblnk(extlaw(j))),j=1,nextlaw)
      else
           write(4,'(11(A,2x))')      "EXTINC_LAW      NONE "
      endif
       if (nmodext .gt. 0) then 
            write(4,'(A,20(I6,1x))') "MOD_EXTINC   ",
     >                                    (modext(i),i=1,nmodext)
       else 
           write(4,'(A,20(I6,1x))') "MOD_EXTINC   ",
     >                                    (modext(i),i=1,2)
       endif 
      write(4,'(A,500(f6.3,1x))')"EB_V             ",(ebv(j),j=1,nebv)
      write(4,'(2A)')     "EM_LINES         ",emlines(1:lnblnk(emlines))
      if (libphys(1:1) .eq. 'Y') then
        write(4,'(A,1x,I6,1x,500(f6.3,1x,2(I12,2x)))') 'Z_RECORD ',iform
     >     ,(zrec(1,k),(idint(zrec(j,k)),j=2,3),k=1,iform) 
        do k = 1,iform
          write(UO,'("Z_RECORD  ",f6.3,i12,1x,i12,2x,i12)')
     >  zrec(1,k),(IDINT(zrec(j,k)),j=2,3),IDINT(zrec(3,k)-zrec(2,k)+1)
        enddo
      endif  
      write(4,'(2A)')     "CREATION_DATE    ",cr_date(1:lnblnk(cr_date))
c  
      close(1)
      close(2)
      close(4)
      if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y") close(5)
      if(UO.eq.30) close(UO)

c 999  format("# EB_V        : ",10(f8.3,1x))
 1010 format(i6,1x,i4,1x,f5.3,1x,E12.6,1x,f6.3,1x,E12.6,1x,E12.6,
     >       1x,i8,1x,i3,1x,1000(E13.6,1x))
 1020 format(i6,1x,i4,1x,f5.3,1x,E12.6,1x,f6.3,1x,E12.6,1x,E12.6,
     >       1x,i8,1x,i3,1x,1500(E13.6,1x))
 1030 format("mod ->",I6," ext law ",i4," extinction ",f6.3,
     >       " age(Gyr) ",E12.6,a,$)
c
      stop
 56   write (UO,*) 'File ',file(1:lnblnk(file)),
     > ' not found -> STOP '
      end
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine colgrid(flux,wave,iw,zoss,fz,
     >     flveg,lbveg,nlveg,lambf,repf,jmax,imag,magtyp,mag)
c
      implicit none
c
      integer*4  i,j,n,iw,imag,wmax,amax,maxsize,nbf
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_wave.decl'
      parameter  (amax=500,maxsize=110000)
      integer*4  nlrep,nlveg,jmax(nbf)
c
      real*8     mag(nbf),h0,om0,l0,zoss,z1,funz,fz
      real*8     rep(wmax),lamb(wmax)
      real*8     wave(wmax),flux(wmax),f(wmax),w(wmax)
      real*8     magab,magveg,abveg,abcor(nbf)
      real*8     lbveg(wmax),flveg(wmax)
      real*8     lambf(nbf,wmax),repf(nbf,wmax)
      character*4096     magtyp

      external  funz
      common /univ/ h0,om0,l0
c
c  put flux and lambda according to z , exti opa...
      z1 = 1+zoss     
c      fz = funz(zoss,h0,om0,l0)
      do n = 1 , iw
         f(n) =  flux(n)/z1                 ! erg.s-1.cm-2.A-1
         if (wave(n).le.50) f(n) = 0.d0
         w(n) = wave(n)*z1 
      enddo
c
      do i =  1, imag
         nlrep=jmax(i)
         do j = 1,nlrep
	   lamb(j) = lambf(i,j)
           rep(j) = repf(i,j)
 	 enddo
c    check if lambda start pass fully through the filter
         if (lamb(nlrep).le.w(iw)) then
c            write(*,*) i,iw,nlrep,nlveg
            call cal_mag(f,w,iw,rep,lamb,nlrep,
     >      flveg,lbveg,nlveg,magtyp,magveg,magab,abveg) 
c            magab=0
c            magveg=0
c            abveg=0 
c            write(*,*) i,magab
            abcor(i) = abveg
            if (magtyp(1:1) .eq. 'A')   mag(i)= magab  + fz 
            if (magtyp(1:1) .eq. 'V')   mag(i)= magveg + fz
          else
            mag(i)   = 99
            abcor(i) = 99
          endif 
      enddo
c
      return 

      end
c




