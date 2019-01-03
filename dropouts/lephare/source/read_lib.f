c
      subroutine read_lib(lib,nlib,size,
     >zstep,zmax,dz,nrec,typlib,
     >modlib,extlib,ebvlib,zlib,agelib,reclib,
     >maglib,klib,emlib,emlines,ldustlib,
     >nzf,zf,zflib,
     >starlib,qsolib,gallib,physpara,
     >extlaw,nextlaw,ebv,nebv,modext,
     >hub,omeg,lbda,magtyp,valf,imagm,filters,nmodel,
     >fobs,minkcol)
c  
c     New :  typlib(size) = 1 for GAL / 2 for QSO / 3 for STAR 
c      program read_lib
      implicit none 
      integer*4     nl,lnblnk,i,j,k,nrec,size,ntest,test,nmodel
      integer*4     nfilt,imagm,vali,reclibmax
      character*4096 zpdir,zpwork,file,valc,typ
c
      integer*4  maxsize
      parameter  (maxsize=110000)
      integer*4     nbf,nused,chisize
      INCLUDE 'out_unit.decl' 
      INCLUDE 'dim_filt.decl'      
      INCLUDE 'dim_zchi.decl' 
      integer*4     nlib,recmax,modext(20),nebv,nzf,nextlaw
      real*8        zstep,zmax,dz,hub,omeg,lbda
      real*8        zstepq,zmaxq,dzq,hubq,omegq,lbdaq
      real*8        ebv(500),zf(500)
      character*4096 lib(5),starlib,gallib,qsolib,filtname(3),filters
      character*4096 magtyp,valf(nbf),extlaw(10)
      character*4096 paravc(500),param,emlines
      integer*4     modlib(size),reclib(size),extlib(size),typlib(size)
      real*8        ebvlib(size),zflib(size)
      real*8        zlib(size),agelib(size),maglib(nbf,size)
      real*8        klib(nbf,size),emlib(nbf,size),ldustlib(size)
c physical parameters 
      integer*4     kmax,kphys,recp,modp
      real*8        physpara(50,maxsize)
c filter choice for Absolute magnitudes  
      real*8        kcolormin(nbf,chisize),kcolormax(nbf,chisize)
      real*8        kcolordiff(nbf,chisize),minkcol(chisize,nbf,nbf)
      real*8        kcolor
      integer*4     l,fobs(chisize,nbf,nbf),fref,indexz,ind

 
c      If you want to give a name to the output screen file
      if(UO.eq.30)
     .      open(UO,file='/tmp/screenReadLib.dat',status='unknown') 
c
cccccccccccccccccccccccccc
c Environemental  Variable 
      call getenv('LEPHAREDIR',zpdir)
      if (lnblnk(zpdir) .eq. 0) then
        write(UO,*) 'WARNING :  variable LEPHAREDIR not defined'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(UO,*) 'WARNING :  variable LEPHAREWORK not defined'
        stop
      endif
c
      write(UO,*) ' nb libraries : ',nlib
cccccccccccccccccccccccccccccccccccccccc
c READING LIBRARIES 
      nrec=0
      ntest=0
      nfilt=0
      zstep=.1
      zmax=1.d0
      dz=0.d0
      nebv=1
      nzf=0
      do k = 1, 500
        ebv(k)= 0.d0
        zf(k) = 0.d0
      enddo  
      typ=' '
      hub=70.d0
      omeg=0.3
      lbda=0.7
      do i = 1, 10 
        extlaw(i)='NONE'
      enddo
      do i = 1, 20 
        modext(i)=0
      enddo
      nmodel=0 
      nused=0
      do nl = 1, nlib   ! loop on libraries used (1->3)
c
c reading doc files for each library 
        file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >            // lib(nl)(1:lnblnk(lib(nl))) // '.doc'
c
        param='LIB_TYPE'
        call read_para2(param,file,paravc,test)
        valc=paravc(1)
c
        param='NUMBER_ROWS'
        call read_para(param,file,paravc,test)
        read(paravc(1),'(i8)') recmax
c
        param='RECORD_LENGTH'
        call read_para(param,file,paravc,test)
        read(paravc(1),'(i10)') reclibmax 
c
        param='MAG_TYPE'
        call read_para2(param,file,paravc,test)
        magtyp=paravc(1)
c
        param='FILTER_FILE'
        call read_para2(param,file,paravc,test)
        filtname(nl)=paravc(1)
        if (nl.gt.1 .and. 
     >    filtname(nl-1)(1:lnblnk(filtname(nl-1))).ne.
     >    filtname(nl)(1:lnblnk(filtname(nl))) ) then 
           write(UO,*) ' Filter file is different '
           write(UO,*) ' between the libraries --> STOP '
        write(UO,*) nl-1,':',filtname(nl-1)(1:lnblnk(filtname(nl-1)))
        write(UO,*) nl,':',filtname(nl)(1:lnblnk(filtname(nl)))
           stop
        endif
        filters=paravc(1)
c
        param='FILTERS'
        call read_para2(param,file,paravc,test)
        imagm=test
        do j = 1,test
           valf(j) = paravc(j)
        enddo
cccccccccccccccccccccccccccccccccccccccccccccc
        if (valc(1:3).eq.'STA') then 
          reclibmax =  8*(3+imagm)
          param='LIB_NAME'
          call read_para2(param,file,paravc,test)
          starlib=paravc(1)
c          write(UO,*) ' read STAR keywords...'
          nused=nused+1
cccccccccccccccccccccccccccccccccccccccccccccc
        elseif (valc(1:3).eq.'GAL') then
c
          param='EM_LINES'
          call read_para(param,file,paravc,test)
          emlines = paravc(1) 
c 
          if (emlines(1:1).eq."Y") then 
              reclibmax =  8*(8+imagm*3)
          else
              reclibmax =  8*(8+imagm*2)
          endif     
c
          param='LIB_NAME'
          call read_para2(param,file,paravc,test)
          gallib = paravc(1) 
c
          param='NUMBER_SED'
          call read_para(param,file,paravc,test)
          read(paravc(1),'(i8)') nmodel        
c
          param='Z_STEP'
          call read_para2(param,file,paravc,test)
          read(paravc(1),'(f8.3)') zstep
          read(paravc(2),'(f8.3)') zmax
          read(paravc(3),'(f8.3)') dz
c
          param='COSMOLOGY'
          call read_para2(param,file,paravc,test)
          read(paravc(1),'(f8.3)') hub
          read(paravc(2),'(f8.3)') omeg
          read(paravc(3),'(f8.3)') lbda
c          write(UO,*) hub,omeg,lbda
c
          param='EXTINC_LAW'
          call read_para2(param,file,paravc,test)
          nextlaw=test
          do i = 1, nextlaw  
             extlaw(i) = paravc(i) 
          enddo
c
          param='EB_V'
          call read_para2(param,file,paravc,test)
          nebv=test
          do i = 1,nebv
             read(paravc(i),'(f8.3)') ebv(i)
          enddo
c
          param='MOD_EXTINC'
          call read_para2(param,file,paravc,test)
          if (test .ge. 2) then 
            do i = 1,test
             read(paravc(i),'(i6)') modext(i)
            enddo
          endif
c
          param='Z_FORM'
          call read_para2(param,file,paravc,test)
          nzf=test
          if (nzf.ge.1 .AND. nzf.eq.nmodel) then
            do i = 1,nzf
               read(paravc(i),'(f8.3)') zf(i)
            enddo
            reclibmax =  8*(6+imagm*2)
          else
             nzf=0
          endif         
c
          nused=nused+2
c
ccccccccccccccccccccccccccccccccccccccccccccccc
        elseif (valc(1:3).eq.'QSO') then
c
          reclibmax =  8*(8+imagm*2)
c
          param='LIB_NAME'
          call read_para2(param,file,paravc,test)
          qsolib = paravc(1) 
c
          param='Z_STEP'
          call read_para2(param,file,paravc,test)
          read(paravc(1),'(f8.3)') zstepq
          read(paravc(2),'(f8.3)') zmaxq
          read(paravc(3),'(f8.3)') dzq
c
          param='COSMOLOGY'
          call read_para2(param,file,paravc,test)
          read(paravc(1),'(f8.3)') hubq
          read(paravc(2),'(f8.3)') omegq
          read(paravc(3),'(f8.3)') lbdaq
c          write(UO,*) hubq,omegq,lbdaq
          nused=nused+4
c
        endif 
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        write(UO,*) ' Library   : ',nl,' ',valc(1:3),
     >' no records:',recmax,reclibmax
c
        if (nl.gt.1 .and. imagm.ne.nfilt) then
           write(UO,*) ' Number of filters is different '
           write(UO,*) ' between the libraries --> STOP '
           stop
        elseif (nl.gt.1 .and. 
     >   typ(1:lnblnk(typ)).ne.magtyp(1:lnblnk(magtyp))) then
           write(UO,*) ' Magnitude type is different '
           write(UO,*) ' between the libraries --> STOP '
           stop   
        endif 
        if (imagm.gt.nbf) then 
          write(UO,*) ' Num filters :',imagm,' default size :',nbf
          write(UO,*) ' -> The dimension of filters exceeds the '
          write(UO,*) ' default s one  in dim_filt.decl  '
          write(UO,*) ' dim_filt.decl must be increased  '
          STOP
        endif   
        ntest=ntest+recmax
        if (ntest.gt.size) then
          write(UO,*) ' dim lib :',ntest,' default size :',size
          write(UO,*) '-> The dimension of library exceeds the'
          write(UO,*) ' default s one  in dim_lib.decl  '
          write(UO,*) ' dim_lib.decl must be increased  '
          STOP
        endif   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c reading the library 
c
        file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     > // lib(nl)(1:lnblnk(lib(nl))) // '.bin'
        open(2,file=file,status='unknown',access='direct',
     >               recl=reclibmax)
        do i = 1, recmax 
          nrec=nrec+1
ccccccccccccccccccccccccccccccccccccccccccccccc
          if (valc(1:3) .eq. "STA") then
c
             typlib(nrec) = 3
c
             read(2,rec=i) modlib(nrec),reclib(nrec),vali,
     >                    (maglib(k,nrec),k=1,vali)
             extlib(nrec) = 0
             ebvlib(nrec) = 0.
             ldustlib(nrec)= -99.
             zlib(nrec)   = 0.
             agelib(nrec) = 0.
             zflib(nrec)  = 0.
             do k = 1,vali 
               emlib(k,nrec)= 0.
               klib(k,nrec) = 0.
             enddo  
cccccccccccccccccccccccccccccccccccccccccccccc
          elseif (valc(1:3) .eq. "QSO") then
c
             typlib(nrec) = 2
c
             read(2,rec=i) modlib(nrec),extlib(nrec),ebvlib(nrec),
     >                ldustlib(nrec),zlib(nrec),agelib(nrec),
     >                reclib(nrec),vali,(maglib(k,nrec),k=1,vali),
     >                     (klib(k,nrec),k=1,vali)
             zflib(nrec)  = 0.
             do k = 1,vali 
               emlib(k,nrec)= 0.
             enddo  
cccccccccccccccccccccccccccccccccccccccccccccc
          elseif (valc(1:3) .eq. "GAL") then
c
            typlib(nrec) = 1
c
            if (nzf.eq.0) then              
              if (emlines(1:1) .eq. "Y") then 
                read(2,rec=i) modlib(nrec),extlib(nrec),ebvlib(nrec),
     >                   ldustlib(nrec),zlib(nrec),agelib(nrec),
     >                   reclib(nrec),vali,(maglib(k,nrec),k=1,vali),
     >               (klib(k,nrec),k=1,vali),(emlib(k,nrec),k=1,vali)
              else
                read(2,rec=i) modlib(nrec),extlib(nrec),ebvlib(nrec),
     >                        ldustlib(nrec),zlib(nrec),agelib(nrec),
     >                   reclib(nrec),vali,(maglib(k,nrec),k=1,vali),
     >                         (klib(k,nrec),k=1,vali)
                do k = 1,vali 
                   emlib(k,nrec)= 0.
                enddo  
              endif
              zflib(nrec) = 0.
            else 
              extlib(nrec)=1
              ldustlib(nrec)=-99.
              do k = 1,vali 
                 emlib(k,nrec)= 0.
              enddo  
              reclib(nrec) = 1
              read(2,rec=i) modlib(nrec),ebvlib(nrec),
     >                    zlib(nrec),agelib(nrec),zflib(nrec),vali,
     >            (maglib(k,nrec),k=1,vali),(klib(k,nrec),k=1,vali)
            endif  
          endif   
c
        enddo
        close(2) 
c
        nfilt=imagm
        typ=magtyp
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   reading physical parameters from gallib in lib_bin
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(valc(1:3).eq.'GAL') then  ! Only for GAL 
  
           file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >          // gallib(1:lnblnk(gallib)) // '.doc'          
           param='NUMBER_ROWS'
           call read_para(param,file,paravc,test)
           read(paravc(1),'(i8)') kmax
c
           file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >         // gallib(1:lnblnk(gallib)) // '.phys'
          write(UO,*) '               and reading Physical parameters'
c     >     file(1:lnblnk(file)),' with',kmax
c    nrec, model, Age, lg(Luv), lg(LR), lg(Lk),  lg(Ltir), Mass, SFR     Met  Tau   D4000
c                 (yr) (  ---- erg/s/Hz ---  )     (Lo)    (Mo)  (Mo/yr)      (yr)  dimless
c              para 1    2       3         4        5       6     7       8     9   10 
           if (kmax .gt. 0 ) then 
              open(2,file=file,status='unknown',err=56)
              k = 0
              do while (.true.)
                 k = k + 1
                read(2,*,end=12) recp,modp,(physpara(i,k),i=1,10)
              enddo  
 12           kphys = k - 1
             close(2)
             if (kphys.ne.kmax) 
     >       write(UO,*) ' Size problem for the  lib_bin library ',
     >                  gallib(1:lnblnk(gallib)),kphys,kmax          
           endif
c           write(UO,*) ' physical parameters read' 
        endif 
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Part for the best filter choice for Mag Abs  
cccccccccccccccccccccccccccccccccccccccccccc
        if(valc(1:3).eq.'GAL') then  ! Only for GAL 
c
          do i = 1 ,indexz(zmax,zstep,dz)
             do j = 1, vali
                do k = 1 , vali 
                   fobs(i,j,k)    = 0
                   minkcol(i,j,k) = 1.d10
                enddo
             enddo
          enddo
c   
          do fref=1,vali                        ! Loop on reference filter
             do i = 1, indexz(zmax,zstep,dz)  
              do k=1, vali
                kcolormin(k,i) =  1.d10
                kcolormax(k,i) = -1.d10
              enddo
             enddo   
c
             do i = 1, recmax                   ! For all library
               do k=1,vali                      ! For each selection filter
c        Find maximum difference in kcolor for each z step and each filter     
                  ind=indexz(zlib(i),zstep,dz)                   !Index in z step
                  kcolor=maglib(fref,i)-maglib(k,i)-klib(fref,i)
                  kcolormin(k,ind)   =  dmin1(kcolor,kcolormin(k,ind))
                  kcolormax(k,ind)   =  dmax1(kcolor,kcolormax(k,ind))
                  kcolordiff(k,ind)  =  dabs(kcolormax(k,ind)
     >                                  -kcolormin(k,ind))
               enddo
             enddo
c
c       Find the adapted selection filter       
             do i=1,indexz(zmax,zstep,dz)       ! Loop on z step
                do j=1,vali                     ! Loop on the prefered order
                   do k=1,vali                  ! Loop on the selection filter
                     do l=j,1,-1  
c                 If this best filter is already found, try another one
                         if(k.eq.fobs(i,fref,l))goto 1 
                     enddo
c        Maximum difference in kcolor
                     if(minkcol(i,fref,j).gt.kcolordiff(k,i))then
                         minkcol(i,fref,j)= kcolordiff(k,i)  
                         fobs(i,fref,j)   = k   ! Best filter ordered in j 
                     endif
 1                   continue
                   enddo
                enddo
             enddo 
c
          enddo !End  loop on reference filter
c
        endif
c
      enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c few tests between Galaxy and QSO
      if ( nused.eq.4 .or. nused.eq. 5 ) then
         zstep=zstepq
         zmax= zmaxq
         dz=   dzq
         hub=  hubq
         omeg= omegq
         lbda= lbdaq
      elseif (nused.ge.6) then   
            
        if (zstepq.ne.zstep.or.zmaxq.ne.zmax.or.dzq.ne.dz.or.
     >     hubq.ne. hub.or.omegq.ne.omeg.or.lbdaq.ne.lbda) then 
           write(UO,*) zstep,zmax,dz,hub,omeg,lbda  
           write(UO,*) zstepq,zmaxq,dzq,hubq,omegq,lbdaq  

           write(UO,*) "Different parameters for z scaling and/or"
           write(UO,*) " cosmological parameters have been used  "
           write(UO,*) " for QSO and GALAXY libraries -> STOP    "
           stop
        endif
      endif  
c
      RETURN
 
 56   write (UO,*) 'File ',file(1:lnblnk(file)),' not found -> STOP '

      END
c

