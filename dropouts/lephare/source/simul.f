c
c     last modif 13/03/01
c     Creation of Random catalogs based on the galaxy libraries  
c     and on an initial luminosity function 
c     The program needs as input : 
c     -1- A magnitude library 
c     -2- An input LF per type and its evolution
c     -3- The behaviour of errors versus magnitude 
c
      PROGRAM  SIMULCAT
      implicit none
c
      integer*4 nbf,inlib,sbmax,nlf,nlfev
      integer*4  nlft
      INCLUDE   'dim_filt.decl'
      INCLUDE   'dim_lib.decl'
      integer   idummy,vali
      parameter (sbmax=310)
      integer*4 i,j,k,l,m,n,p,lref,nobj
      integer*4 kmax,nmagmax,nzmax,nztrans,lf,imagm,bref                 
      integer*4 nfilt,iref,rnum,lfref,lnblnk
      integer*4 linemod(2,inlib),model(inlib)
      integer*4 lftyp(3,100),tottyp(100),lfevt(100),lfprof(100)
      integer*4 nexti,lexti(100)
      real*8    mag(nbf,inlib),z(inlib),kcor(nbf,inlib)
      real*8    zfor(inlib),exti(inlib),age(inlib)
      real*8    ageref,ebvref,zforref
      real*8    magref(nbf),magrefo(nbf),magref1(nbf)
      real*8    flran(nbf),errflran(nbf),error_flux,bidf(nbf)
      real*8    uplim1(nbf),uplim(nbf)
      real*8    errmag(nbf),mabsz(nbf)
      real*8    Ae(nbf),Be(nbf),Ce(nbf),De(nbf)
      real*8    flmoy(nbf),flwidth(nbf)
      real*8    scale,color,temp
      real*8    maglim,sig,sigap,dmag,mabs,mags,magso
      real*8    val,error,sig_noise
      real*8    x,dl,c,h0,om0,oml,dz,dzsup,zbin,zmax,zmod
      real*8    zref,kcref,kciref
      real*8    exp_dens
      real*8    surf,nexpd,nexp
      real*8    hzsteidel,hz3lf,ecor,phicor,alpcor,zeta,etal,etap,etas
      real*8    lftypz(2,100),lftypm(2,100)
      real*8    lfevz(100),lfphi(100),lfmst(100),lfalp(100)
      real*8    phi0(100),mst0(100),alpha0(100)
      real*8    fcorr(nbf)
      integer*4 ntyp
c SB effect parameters
      integer*4 is,ib,ismax,ibmax,val1,valt,typ,nsig
      integer*4  pm,pmp,nbmag,compmag(10),nbaper,sigtmax(2)
      integer*4 imin,imax
      real*8    mulim(nbf)
      real*8    val2,val3,val4,val5
      real*8    sigt(2,sbmax,sbmax),beta(2,sbmax,sbmax)
      real*8    gti(2,sbmax,sbmax),lgt(2,sbmax,sbmax)
      real*8    sigtr(sbmax),betar(sbmax),lgtr(sbmax),gtir(sbmax)
      real*8    re0,re_mabs,sigmat,da,fac,pi,mabs_bv,mabs_bd,magp
      real*8    dmin,mag_fac,ke,sblim,fdmin(nbf)
      real*8    rad_scale(10),eta_p,rad_aper(5)
      real*8    seeing(nbf),re0as
      real*8    magsb(10,nbf),errmagsb(10,nbf)
      real*8    thetal,theta(nbf),theta_mag
      real*8    thetamag(10,nbf)
      real*8    colz0(nbf),colorz0(nbf,inlib),brvega,brcol(inlib)
c
c ERROR SIM parameters 
      real*8    mag_ran,sn,zp(nbf),gain(nbf),pixscale(nbf),abcor(nbf)
      real*8    fuge_fac(nbf),mag_lim,uplima(10,nbf),uplima1(10,nbf)
c
c other parameters 
      real       poidev
      real*4     ran1,nexpo
      real*8     schechter
      real*8     min,max,zmin
c
      character*512  inf,simf,errf,lff,lfevol,convf
      character*512  type,sb,sexerr
      character*512  cr_date,fdate,mag_typ
c
      real*8        paravr(100)  
      integer*4     paravi(100)
      character*512 paravc(100),namefl(nbf),nameff(nbf)
      character*1024 str

c
      character*512 zpdir,zpwork,config,param,filters
      character*512 outpara,wpara(100),str_out(100)
      integer*4 iwpara,iwout
c
      integer*4  test,jmax
c
      common /dum/        idummy
      common /errorfile/  errf
      EXTERNAL schechter,sig_noise,maglim,error_flux
      external x,exp_dens,error,poidev,hzsteidel,hz3lf,re_mabs,ran1
cccccccccccccccccccccccccccccccccc
c initialisation 
      lref=0
      bref=0
      color=0
      nbmag=0
      nbaper=0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c constant parameters
      c=300000.
      pi=3.141592654
      fac=180*3600/pi
ccccccccccccccccccccccccccccccccccccc
c  Environemental variable 
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
c
ccccccccccccccccccccccccccccccccccccc
c  Input parameter from config file
      param='-c'
      call get_conf(param,config,test)
      if (test.ne.1) call err_option(param,1)
      call get_path(config)
c      write(*,'(A)') config(1:lnblnk(config)) 

      write(6,'(2(A,1x))') '   >>> reading keywords from ',
     > config(1:lnblnk(config))
c
cccccccccccccccccccccccccccccccccccccc
c  Reading Filter FILE 
      param='-FILTER_FILE'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  filters=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
      write(6,'(A)') '   >>> filter characteristics ...'
      call zeropoint(filters,zp,abcor,flmoy,flwidth,fcorr,nfilt)
c  Extract names of filters 
        filters = zpwork(1:lnblnk(zpwork)) //'/filt/'
     >            // filters(1:lnblnk(filters))
        open(1,file=filters,status='old')
        read(1,'(A)') str
        call val_string(str,paravc,test)
        read(paravc(2),'(i10)') nfilt
        do i =  1, nfilt
          read(1,'(A)') str
          call val_string(str,paravc,test)
          read(paravc(2),'(i10)') jmax
          nameff(i) = paravc(3)
          do j = 1,jmax
             read(1,'(A)') str
          enddo   
       enddo
       CLOSE(1)
cccccccccccccccccccccccccccccccccccccccccccccc
c
      param='-GAL_LIB_OUT'
      call getc_option(param,config,100,paravc,test)
      if (test.eq.1 ) then
           inf=paravc(1)(1:lnblnk(paravc(1)))
      else
         call err_option(param,1)         
      endif      
c
      param='-LF_FILE'
      call getc_option(param,config,100,paravc,test)
      if (test.eq.1 ) then
           lff=paravc(1)(1:lnblnk(paravc(1)))
           call get_path(lff)
      else
         call err_option(param,1)         
      endif      
c
      param='-LF_EVOL'
      call getc_option(param,config,100,paravc,test)
      if (test.eq.1 ) then
           lfevol=paravc(1)(1:lnblnk(paravc(1)))
           call get_path(lfevol)
      else
         lfevol='NONE'
         call err_option(param,2)         
      endif    
c      
      param='-ERROR_FILE'
      call getc_option(param,config,100,paravc,test)
      if (test.eq.1 ) then
           errf=paravc(1)(1:lnblnk(paravc(1)))
           call get_path(errf)
      else
         call err_option(param,1)         
      endif      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Read error behaviours 
      write(*,*) '   >>> reading ',errf(1:lnblnk(errf))
      open(1,file=errf(1:lnblnk(errf)),status='old') 
      read(1,*) 
      i=0
      do while (.true.)
         i=i+1
         read(1,*,end=12) val,Ae(i),Be(i),Ce(i),De(i)
      enddo
 12   close(1)
      nlf = i -1 
      if (nlf.ne.nfilt) then
         write(*,*) ' Number of filters in ',errf(1:lnblnk(errf))
         write(*,*) ' NOT EQUAL TO ',filters(1:lnblnk(filters))
         STOP
      endif   
c
      param='-OUTPUT_FILE'
      call getc_option(param,config,100,paravc,test)
      if (test.eq.1 ) then
           simf=paravc(1)(1:lnblnk(paravc(1)))
           call get_path(simf)
      else
         call err_option(param,1)         
      endif   
c
      write(6,*) '   >>> reading output parameter file ... '
      param='-PARA_OUT'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1) outpara=paravc(1)(1:lnblnk(paravc(1)))
      if (test.eq.1) call get_path(outpara)
      if (test.ne.1) call err_option(param,1)
c  reading output parameter
      open(1,file=outpara,status='unknown')
      iwpara = 0 
      do while (.true.)
         read(1,'(A)',end=14) str 
         if (str(1:1) .ne. '#' .AND. str(1:1) .ne. ' ' 
     >  .AND. str(1:1) .ne. char(9) .AND. str(1:1) .ne. char(13)) then
           call val_string(str,paravc,test)
           iwpara = iwpara + 1
           wpara(iwpara) = paravc(1)
c           write(6,*) iwpara,wpara(iwpara)(1:lnblnk(wpara(iwpara)))
         endif
      enddo
 14   close(1)
c   
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Read parameters 
      param='-M_STAR'
      call getf_option(param,config,100,paravr,test)
      nlft=test
      if (test.gt.0) then
         do i = 1,nlft
            mst0(i)=paravr(i)
         enddo
      else
         call err_option(param,3)  
      endif  

      param='-PHI_STAR'
      call getf_option(param,config,100,paravr,test)
      if (test.ne.nlft) then
          call err_option(param,3)
      else
         do i = 1,nlft
            phi0(i)=paravr(i)
         ENDDO
      endif  

      param='-ALPHA'
      call getf_option(param,config,100,paravr,test)
      if (test.ne.nlft) then
          call err_option(param,3)
      else
         do i = 1,nlft
            alpha0(i)=paravr(i)
         enddo
      endif  
c
      param='-LF_PROFILE'
      call geti_option(param,config,100,paravi,test)
      if (test.ne.nlft) then
          call err_option(param,3)
      else
         do i = 1,nlft
            lfprof(i)=paravi(i)
         enddo
      endif        
c
      param='-ETAL_EVOL'
      call getf_option(param,config,100,paravr,test)
      if (test.eq.1) then 
         etal=paravr(1)
      else
         etal=0.
      endif
      param='-ETAP_EVOL'
      call getf_option(param,config,100,paravr,test)
      if (test.eq.1) then 
         etap=paravr(1)
      else
         etap=0.
      endif
      param='-ETAS_EVOL'
      call getf_option(param,config,100,paravr,test)
      if (test.eq.1) then 
         etas=paravr(1)
      else
         etas=0.
      endif

c test if LF EVOL file exist -> cancels eta 
      if ( lfevol(1:4) .ne. 'NONE')  then 
         etal = 0 
         etap = 0 
         etas = 0 
      endif
c
      param='-ZETA_EVOL'
      call getf_option(param,config,100,paravr,test)
      if (test.eq.1) then 
         zeta=paravr(1)
      else
         zeta=0.
      endif
c 
      param='-MAG_LIMIT'
      call getf_option(param,config,100,paravr,test)
      
      if (test.eq.2) then
         min=paravr(1)
         max=paravr(2)
      else
         call err_option(param,3)  
      endif         
c
      param='-Z_LIMIT'
      call getf_option(param,config,100,paravr,test)
      if (test.eq.1) then
         zmax=paravr(1)
      else
         call err_option(param,3)  
      endif         
c
      param='-SIGMA_UPLIM'
      call getf_option(param,config,100,paravr,test)      
      if (test.eq.1) then
         sig=paravr(1)
      else
         call err_option(param,3)  
      endif         
c
      param='-AREA_SIM'
      call getf_option(param,config,100,paravr,test)
      
      if (test.eq.1) then
         surf=paravr(1)
      else
         call err_option(param,3)  
      endif         
c
      param='-MAG_REF'
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  iref=paravi(1)
      if (test.ne.1)  call err_option(param,3)
c
      param='-LF_REF'
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  lfref=paravi(1)
      if (test.ne.1)  call err_option(param,3)
c
      param='-DUMMY_INT'
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  idummy=paravi(1)
      if (test.ne.1)  idummy=-54396523
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c  PARAMETERS FOR SURFACE BRIGHTNESS EFFECTS cccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
      param='-SB_EFFECT'
      call getc_option(param,config,100,paravc,test)
      if (test.eq.1 ) then
           sb=paravc(1)(1:lnblnk(paravc(1)))
      else
         call err_option(param,1)         
      endif      
      if (sb(1:1).eq.'Y' .or. sb(1:1).eq.'y' ) then 
c
         param='-CONV_FILE'
         call getc_option(param,config,100,paravc,test)
         if (test.eq.1 ) then
            convf=paravc(1)(1:lnblnk(paravc(1)))
            call get_path(convf)
         else
            call err_option(param,1)         
         endif      
c
         param='-FILT_B'
         call geti_option(param,config,1,paravi,test)
         if (test.eq.1)  bref=paravi(1)
         if (test.ne.1)  call err_option(param,3)
c
         nbmag = 1       ! at least isophotal 
         compmag(1) = 1  ! at least isophotal
c
         param='-PSEUDO_RAD'
         call getf_option(param,config,100,paravr,test)      
         if (test.eq.1) then
            rad_scale(1)=paravr(1)
            nbmag = nbmag + test
            compmag(2) = 1
         else
            compmag(2)  = 0
            rad_scale(1)= 0
            call err_option(param,2)  
         endif         
c
         param='-AUTO_RAD'
         call getf_option(param,config,100,paravr,test)      
         if (test.eq.1) then
            rad_scale(2)=paravr(1)
            nbmag = nbmag + test 
            compmag(3) = 1
         else
            compmag(3)  = 0
            rad_scale(2)= 0
            call err_option(param,2)  
         endif         
c
         param='-ETA_PETRO'
         call getf_option(param,config,100,paravr,test)      
         if (test.eq.1) then
            eta_p=paravr(1)
            nbmag = nbmag + test
            compmag(4) = 1
         else
            compmag(4) = 0
            call err_option(param,2)  
         endif         
c
         param='-APER_DIAM'
         call getf_option(param,config,100,paravr,test)  
         nbaper=test
         if (test.ge.1) then
            do i = 1, test 
               rad_aper(i)=paravr(i)
               compmag(4+i) = 1 
            enddo  
            nbmag = nbmag + test
         else
            do i = 5,10
               compmag(i) = 0 
            enddo   
            call err_option(param,2)  
         endif         
         write(*,*) (compmag(i),i=1,9)
c
         param='-SEEING'
         call getf_option(param,config,100,paravr,test)
         nlf=test
         if (test.ge.1) then
           do i = 1,nlf
             seeing(i) =  paravr(i)
           enddo
         else
            call err_option(param,1) 
         endif            
         if (nlf.ne.nfilt) then
           write(*,*) ' Number of filters in SEEING'
           write(*,*) ' NOT EQUAL TO ',filters(1:lnblnk(filters))
           STOP
         endif   
c
         param='-DMIN'
         call getf_option(param,config,100,paravr,test)
         nlf=test      
         if (test.ge.1) then
            do i = 1,nlf
               fdmin(i)=paravr(i)
            enddo
         else
            call err_option(param,1)  
         endif         
         if (nlf.ne.nfilt) then
           write(*,*) ' Number of filters in DMIN'
           write(*,*) ' NOT EQUAL TO ',filters(1:lnblnk(filters))
           STOP
         endif            
c
         param='-SB_LIM'
         call getf_option(param,config,100,paravr,test)
         nlf=test
         if (test.ge.1) then
           do i = 1,nlf
             mulim(i) =  paravr(i)
           enddo
         else
            call err_option(param,1) 
         endif
         if (nlf.ne.nfilt) then
           write(*,*) ' Number of values in SB_LIM'
           write(*,*) ' NOT EQUAL TO ',inf(1:lnblnk(inf))
           STOP
        endif   
c  TYPE OF ERRORS        
        param='-SEX_ERROR'
        call getc_option(param,config,100,paravc,test)
        if (test.eq.1 ) then
           sexerr=paravc(1)(1:lnblnk(paravc(1)))
        else
           call err_option(param,1)         
        endif      
        if (sexerr(1:1).eq.'Y' .or. sexerr(1:1).eq.'y') then 
c
           param='-PIXEL_SCALE'
           call getf_option(param,config,100,paravr,test)
           nlf=test
           if (test.ge.1) then
             do i = 1,nlf
               pixscale(i) = paravr(i)
             enddo
           else
             call err_option(param,2) 
           endif     
           if (nlf.ne.nfilt) then
             write(*,*) ' Number of values in PIXEL_SCALE'
             write(*,*) ' NOT EQUAL TO ',inf(1:lnblnk(inf))
             STOP
           endif   
c
           param='-ZP_OBS'
           call getf_option(param,config,100,paravr,test)
           nlf=test
           if (test.ge.1) then
             do i = 1,nlf
               zp(i) =  paravr(i)
             enddo
           else
             call err_option(param,2) 
           endif   
           if (nlf.ne.nfilt) then
             write(*,*) ' Number of values in ZP_OBS'
             write(*,*) ' NOT EQUAL TO ',inf(1:lnblnk(inf))
             STOP
           endif   
c
           param='-GAIN'
           call getf_option(param,config,100,paravr,test)
           nlf=test
           if (test.ge.1) then
             do i = 1,nlf
               gain(i) = paravr(i)
             enddo
           else
             call err_option(param,2) 
           endif     
           if (nlf.ne.nfilt) then
             write(*,*) ' Number of values in GAIN'
             write(*,*) ' NOT EQUAL TO ',inf(1:lnblnk(inf))
             STOP
           endif   
c
           param='-FUGE_FACTOR'
           call getf_option(param,config,100,paravr,test)
           nlf=test
           if (test.ge.1) then
             do i = 1,nlf
               fuge_fac(i) =  paravr(i)
             enddo
           else
             call err_option(param,2) 
           endif
           if (nlf.ne.nfilt) then
             write(*,*) ' Number of values in FUGE_FACTOR'
             write(*,*) ' NOT EQUAL TO ',inf(1:lnblnk(inf))
             STOP
           endif
        endif     ! sexerr :  SEXtractor  errors    
      endif       ! sb     :  Surface Brightness effects 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  READING THE PHOTOMETRIC LIBRARY and put magnitudes in AB for LF  
c
      write(6,'(A)') '    >>> reading galaxy library  ...'
      call read_libsim(inf,inlib,dz,zmod,dzsup,h0,om0,oml,
     >                 model,z,mag,kcor,zfor,exti,age,
     >                 namefl,type,imagm,kmax)
c test number of filters 
      if (nfilt.ne.imagm) then
         write(6,*) 'WARNING : Number of filters is not equal between'
         write(6,*) filters(1:lnblnk(filters)),' and '
     >    ,inf(1:(lnblnk(inf)))
         write(6,*) ' check in ',config(1:lnblnk(config)),
     >                 ' the file name in FILTER_FILE -> STOP'
         stop
      endif 
c test name of filters 
      do i = 1, imagm 
         if (nameff(i)(1:lnblnk(nameff(i))).ne.
     >       namefl(i)(1:lnblnk(namefl(i))) ) then
           write(6,*) ' WARNING : Name of filters is not equal '
           write(6,*) ' between filter file and library -> STOP'
           stop
         endif
      enddo   
c
c      zmod=0
      brvega=0
c
      do k = 1, kmax        ! info from input library
c  convert magnitudes from library in AB if needed  
         do j = 1 , nfilt
            if (mag(j,k).ge.99) then 
               mag(j,k)=99.
            else
               if (type(1:1).eq.'V' .or.type(1:1).eq.'v') then 
                  mag(j,k)= mag(j,k) + abcor(j) 
               endif   
            endif
            if (kcor(j,k).ge.99) then 
               kcor(j,k)=99.
            endif
         enddo
c  index range for each model
         if (k.eq.1) linemod(1,model(k)) = k
         if (k.gt.1.and.model(k).gt.model(k-1)) then
            linemod(2,model(k-1))= k-1
            linemod(1,model(k))  = k
         endif   
c 
c         if (z(k).gt.zmod) zmod=z(k)
c  compute the color at z=0 for each galaxy SED with respect to lf_ref band 
         if (z(k).le.1.e-5) then
c  if Surface Brightness effect :
c  vega color@z=0 (B-iref)  for  the SB analysis (re vs Mabs(B))
            if (sb(1:1).eq.'y' .or. sb(1:1).eq.'Y') then
              brvega=mag(bref,k)-abcor(bref)-(mag(iref,k)-abcor(iref))
            endif
c  ab color@z=0: (Band-iref) in AB for all filters
            do j = 1, nfilt
               colz0(j)= mag(j,k) - mag(iref,k)
            enddo  
         endif   
         if (sb(1:1).eq.'y' .or. sb(1:1).eq.'Y') brcol(k)=brvega
         do j = 1, nfilt
            colorz0(j,k)=colz0(j)
         enddo   
      enddo
c
      linemod(2,model(kmax))=kmax
      if (zmod.le.zmax) zmax=zmod
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  read Type for each LF 
      write(*,*) '   >>> reading ',lff(1:lnblnk(lff))
      open(1,file=lff(1:lnblnk(lff)),status='old')
      read(1,*)  
      i=0
      do while (.true.)
         i = i + 1
         read(1,*,end=22) (lftyp(j,i),j=1,3),(lftypz(j,i),j=1,2)
     >                   ,(lftypm(j,i),j=1,2)
         tottyp(i) = lftyp(3,i)-lftyp(2,i) + 1
      enddo   
 22   close(1)
      nlf = i - 1     ! Number of LF types 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  read LF parameters and evolution 
      if ( lfevol(1:4) .eq. 'NONE' ) then
         lfevt(1)=1
         lfevz(1)=0
         lfphi(1)=0
         lfmst(1)=0
         lfalp(1)=0
         nlfev=0
      else   
         write(*,*) '   >>> reading ',lfevol(1:lnblnk(lfevol))
         open(1,file=lfevol(1:lnblnk(lfevol)),status='old')
         read(1,*)  
         i=0
         do while (.true.)
           i = i + 1
           read(1,*,end=26) lfevt(i),lfevz(i),
     >     lfphi(i),lfmst(i),lfalp(i)
           write(*,*) lfevt(i),lfevz(i),
     >     lfphi(i),lfmst(i),lfalp(i)
         enddo   
 26      close(1)
         nlfev = i - 1
         if (nlfev.gt.100) then
            write(6,'(I4,A)') nlfev,
     >                  ' : SIZE greater than the declaration for ',
     >                    lfevol(1:lnblnk(lfevol))
            stop
         endif
      endif   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  surface brigthness included 
      if (sb(1:1).eq.'y' .or. sb(1:1).eq.'Y') then 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c read profile_psf
        write(*,*) '   >>> reading ',convf(1:lnblnk(convf))
        open(1,file=convf(1:lnblnk(convf)),status='old')
        i=0
        is=0
        ib=0
        val=-1
        valt=-1
        sigtmax(1)=0
        sigtmax(2)=0
        do while (.true.)
          i = i + 1
          read(1,*,end=2) val1,val2,val3,val4,val5
          if (val1.ne.valt) then
            is=0
            ib=0
            valt=val1
            val=-1
          endif   
          if (val3.gt.val) then
            val=val3
            is=is+1
            ib=0
          endif   
          ib=ib+1
          if (val1.eq.1) then
            sigt(1,is,ib) = val3
            beta(1,is,ib) = val2
            lgt(1,is,ib)  = val4
            gti(1,is,ib)  = val5
            if (val3.gt.sigtmax(1)) sigtmax(1)=val3
          else
            sigt(2,is,ib) = val3
            beta(2,is,ib) = val2
            lgt(2,is,ib)  = val4
            gti(2,is,ib)  = val5
            if (val3.gt.sigtmax(2)) sigtmax(2)=val3
          endif   
        enddo
 2      close(1)
        ismax=is
        ibmax=ib
        write(*,*) 'i,is,ib,stmax(2):',i,is,ib,sigtmax(1),sigtmax(2)
      endif  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  magnitude range 
      dmag=0.10
      nmagmax=DINT((max-min)/dmag)+1
c  redshift range 
      zmin=0.
      if (zmax.le.6) then
        nzmax  = IDNINT(zmax/dz) + 1
        nztrans= nzmax
      else
        nztrans= IDNINT(6./dz) + 1 
        nzmax  = nztrans + IDNINT((zmax-6.)/dzsup)
      endif  
c  computes the UPPER-LIMITS in 2*FWHM  at  n sigma for std magnitudes
      write(*,*) '   >>> computing upper-limits ...'
      do k = 1, nfilt
         mulim(k) = mulim(k) -abcor(k)  ! put mulim in VEGA for MAGLIM_SIM
         if (sexerr(1:1).eq.'Y' .or. sexerr(1:1).eq.'y') then  
            thetal=2*seeing(k)
            call MAGLIM_SIM(thetal,sig,k,mulim,abcor,zp,
     >                       gain,pixscale,fuge_fac,mag_lim)
            uplim(k) = mag_lim
            call MAGLIM_SIM(thetal,1.d0,k,mulim,abcor,zp,
     >                       gain,pixscale,fuge_fac,mag_lim)
            uplim1(k) = mag_lim

         else   
            uplim(k) = maglim(k,sig,Ae,Be,Ce,De)
            uplim1(k) = maglim(k,1.d0,Ae,Be,Ce,De)
         endif   
c  For apertures : used sig=1
         sigap = 1. 
         do i = 1, nbaper
            thetal=rad_aper(i)
           if (sexerr(1:1).eq.'Y' .or. sexerr(1:1).eq.'y') then  
              call MAGLIM_SIM(thetal,sigap,k,mulim,abcor,zp,
     >                       gain,pixscale,fuge_fac,mag_lim)
              uplima(i,k) = mag_lim

              call MAGLIM_SIM(thetal,1.d0,k,mulim,abcor,zp,
     >                       gain,pixscale,fuge_fac,mag_lim)
              uplima1(i,k) = mag_lim

           else   
              uplima(i,k) = maglim(k,sigap,Ae,Be,Ce,De)  
              uplima1(i,k)= maglim(k,1.d0,Ae,Be,Ce,De)  
           endif 
         enddo   
         mulim(k) = mulim(k) +abcor(k) ! put back mulim in AB 
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  writing input parameters on screen 
       write(6,'(A)')"##############################################"
       write(6,'(A)')"#          Building Mock catalog             #"
       write(6,'(A)')"# with ESO-Sculptor LF in Rvega (Omega_m=1)  #"
       write(6,'(A)')"#     and Magnitudes in AB system            #"
       write(6,'(2A)') "#  LIBRARY FILE         : ",inf(1:lnblnk(inf))
       write(6,'(20(A,1x))')"#  FILTER NAMES         : ",
     >                    (nameff(k)(1:lnblnk(nameff(k))),k=1,nfilt)
       write(6,'(A,1x,100(f6.3,1x))')
     >               "#  Vega-to-AB system    :",(abcor(k),k=1,nfilt)
       write(6,'(A,1x,3(f7.2,1x))')
     >                 "#  Cosmology            : ",h0,om0,oml
       write(6,'(2A)') "#  LF TYPE FILE         : ",
     >                                     lff(1:lnblnk(lff))
       write(6,'(A,10(f8.3,1x))')
     >                 "#  M*   (@z=0)          : ",(mst0(i),i=1,nlft)
       write(6,'(A,10(E12.6,1x))')
     >                 "#  Phi* (@z=0)          : ",(phi0(i),i=1,nlft)
       write(6,'(A,10(f6.2,1x))')
     >                 "#  Alpha(@z=0)          : ",(alpha0(i),i=1,nlft)
       write(6,'(A,100(i4,1x))')
     >                  "# LF_PROFILE           : ",(lfprof(i),i=1,nlft)       
       write(6,'(2A)') "#  LF EVOLUTION FILE    : ",
     >                                     lfevol(1:lnblnk(lfevol))
       write(6,'(A,1x,4(f6.3,1x))') 
     >         "# ETAL-P-S, ZETA PARAMETERS  : ",etal,etap,etas,zeta
       write(6,'(2A)') "#  ERROR FILE           : ",
     >                                     errf(1:lnblnk(errf))
       write(6,'(A,1x,I6)')"#  LF REFERENCE  BAND   : ",lfref
       write(6,'(A,1x,I6)')"#  REFERENCE BAND       : ",iref
       write(6,'(A,1x,3(f6.3,1x))')
     >               "#  Magnitude range, dmag: ",min,max,dmag
       write(6,'(A,1x,4(f6.3,1x))')
     >               "#  Redshift range, dz   : ",zmin,zmax,dz,dzsup
       write(6,'(2A)') 
     >               "#  SURFACE BRIGHTNESS   : ",sb(1:3)
       if (sb(1:1).eq."y" .or. sb(1:1).eq."Y") then
          write(6,'(A,1x,50(f6.2,1x))')
     >               "#  Surf. Bright.(AB)    : ",(mulim(i),i=1,nfilt)
          write(6,'(A,1x,50(f6.2,1x))')
     >               "#  seeing   (arcsec)    : ",(seeing(i),i=1,nfilt)        
          write(6,'(A,1x,50(f6.2,1x))')
     >               "#  Dmin(frac. of seeing): ",(fdmin(i),i=1,nfilt)  
          write(6,'(2A)')
     >               "#  SEXTRACTOR ERRORS    : ",
     >                   sexerr(1:lnblnk(sexerr))  
          if (sexerr(1:1).eq.'Y' .or. sexerr(1:1).eq.'y') then  
             write(6,'(A,1x,50(f6.3,1x))')
     >          "#  sigma, up.lim.(2FWHM):",sig,(uplim(k),k=1,nfilt)
          else
             write(6,'(A,1x,50(f6.3,1x))')
     >          "#  sigma, up. lim.      :",sig,(uplim(k),k=1,nfilt)
          endif  
          if (compmag(1).eq.1)   write(6,'(A,1x,f6.3)')
     >               "#  Mag ISOPHOTAL  with  : ", 1.0 
          if (compmag(2).eq.1)   write(6,'(A,1x,f6.3)')
     >               "#  Mag PSEUDO-ISO with  : ", rad_scale(1) 
          if (compmag(3).eq.1)   write(6,'(A,1x,f6.3)')
     >               "#  Mag AUTOMATIC  with  : ", rad_scale(2) 
          if (compmag(4).eq.1)   write(6,'(A,1x,f6.3)')
     >               "#  Mag PETROSIAN    eta : ", eta_p 
          if (compmag(5).eq.1)   write(6,'(A,1x,5(f6.3,1x))')
     >               "#  Mag APERTURE   with  : ",
     >                  (rad_aper(i),i=1,nbaper) 
          do i = 1, nbaper
            write(6,'(A,1x,50(f6.3,1x))')
     >       "#  sigma, up.lim. Aper  :",sigap,(uplima(i,k),k=1,nfilt)
          enddo
       else
          write(6,'(A,1x,50(f6.3,1x))')
     >          "#  sigma, up. lim.      :",sig,(uplim(k),k=1,nfilt)         
       endif   
       write(6,'(2A)') 
     >               "#  OUTPUT FILE          : ",simf(1:lnblnk(simf))  
       write(6,'(2A)') "# OUTPUT PARAMETER FILE : ",
     >                                     outpara(1:lnblnk(outpara))
       write(6,'(A,1x,f9.6)')
     >               "#  total area (deg2)    : ",surf
       write(6,'(A)')"#                                            #"
       write(6,'(A)')"##############################################"
c
      open(1,file=simf) 
c writing parameter  in output catalog 
       cr_date=fdate()
       write(1,'(A)')"##############################################"
       write(1,'(A)')"#          Building Mock catalog             #"
       write(1,'(A)')"# with ESO-Sculptor LF in Rvega (Omega_m=1)  #"
       write(1,'(A)')"#     and Magnitudes in AB system            #"
       write(1,'(2A)') "#  LIBRARY FILE         : ",inf(1:lnblnk(inf))
       write(1,'(20(A,1x))')"#  FILTER NAMES         : ",
     >                    (nameff(k)(1:lnblnk(nameff(k))),k=1,nfilt)
       write(1,'(A,1x,100(f6.3,1x))')
     >                  "#  Vega-to-AB system    :",(abcor(k),k=1,nfilt)
       write(1,'(A,1x,3(f7.2,1x))')
     >                  "#  Cosmology            : ",h0,om0,oml
        write(1,'(2A)') "#  LF TYPE FILE         : ",
     >                                     lff(1:lnblnk(lff))
       write(1,'(A,50(f8.3,1x))')
     >                  "#  M*   (@z=0)          : ",(mst0(i),i=1,nlft)
       write(1,'(A,50(E12.6,1x))')
     >                  "#  Phi* (@z=0)          : ",(phi0(i),i=1,nlft)
       write(1,'(A,50(f6.2,1x))')
     >                  "#  Alpha(@z=0)          : ",
     >               (alpha0(i),i=1,nlft)
      write(1,'(A,50(i4,1x))')
     >                  "#  LF_PROFILE           : ",
     >               (lfprof(i),i=1,nlft)

       write(1,'(2A)')  "#  LF EVOLUTION FILE    : ",
     >                                     lfevol(1:lnblnk(lfevol))
       write(1,'(A,1x,4(f6.3,1x))') 
     >         "# ETAL-P-S, ZETA PARAMETERS  : ",etal,etap,etas,zeta
       write(1,'(2A)')  "#  ERROR FILE           : ",
     >                                     errf(1:lnblnk(errf))
       write(1,'(A,1x,I6)')"#  LF REFERENCE  BAND   : ",lfref
       write(1,'(A,1x,I6)')"#  REFERENCE BAND       : ",iref
       write(1,'(A,1x,3(f6.3,1x))')
     >               "#  Magnitude range, dmag: ",min,max,dmag
       write(1,'(A,1x,4(f6.3,1x))')
     >               "#  Redshift range, dz   : ",zmin,zmax,dz,dzsup
       write(1,'(2A)') 
     >               "#  SURFACE BRIGHTNESS   : ",sb(1:3)
       if (sb(1:1).eq."y" .or. sb(1:1).eq."Y") then
          write(1,'(A,1x,50(f6.2,1x))')
     >               "#  Surf. Bright.(AB)    : ",(mulim(i),i=1,nfilt)
          write(1,'(A,1x,50(f6.2,1x))')
     >               "#  seeing   (arcsec)    : ",(seeing(i),i=1,nfilt)        
          write(1,'(A,1x,50(f6.2,1x))')
     >               "#  Dmin(frac. of seeing): ",(fdmin(i),i=1,nfilt) 
          write(1,'(2A)')
     >               "#  SEXTRACTOR ERRORS    : ",
     >                   sexerr(1:lnblnk(sexerr))  
          if (sexerr(1:1).eq.'Y' .or. sexerr(1:1).eq.'y') then  
             write(1,'(A,1x,50(f6.3,1x))')
     >          "#  sigma, up.lim.(2FWHM):",sig,(uplim(k),k=1,nfilt)
          else
             write(1,'(A,1x,50(f6.3,1x))')
     >          "#  sigma, up. lim.      :",sig,(uplim(k),k=1,nfilt)
          endif  
          if (compmag(1).eq.1)   write(1,'(A,1x,f6.3)')
     >               "#  Mag ISOPHOTAL  with  : ", 1.0 
          if (compmag(2).eq.1)   write(1,'(A,1x,f6.3)')
     >               "#  Mag PSEUDO-ISO with  : ", rad_scale(1) 
          if (compmag(3).eq.1)   write(1,'(A,1x,f6.3)')
     >               "#  Mag AUTOMATIC  with  : ", rad_scale(2) 
          if (compmag(4).eq.1)   write(1,'(A,1x,f6.3)')
     >               "#  Mag PETROSIAN  with  : ", eta_p 
          if (compmag(5).eq.1)   write(1,'(A,1x,5(f6.3,1x))')
     >               "#  Mag APERTURE   with  : ",
     >                  (rad_aper(i),i=1,nbaper) 
           do i = 1, nbaper
            write(1,'(A,1x,50(f6.3,1x))')
     >       "#  sigma, up.lim. Aper  :",sigap,(uplima(i,k),k=1,nfilt)
          enddo
             else
          write(1,'(A,1x,50(f6.3,1x))')
     >          "#  sigma, up. lim.      :",sig,(uplim(k),k=1,nfilt)         
       endif   
       write(1,'(2A)') 
     >               "#  OUTPUT FILE          : ",simf(1:lnblnk(simf))  
       write(1,'(2A)') "# OUTPUT PARAMETER FILE : ",
     >                                     outpara(1:lnblnk(outpara))
       write(1,'(A,1x,f9.6)')
     >               "#  total area           : ",surf
       write(1,'(2A)')
     >         "#  date of creation     : ",cr_date(1:lnblnk(cr_date)) 
       write(1,'(A)')
     >         "#  Output Format        : "
      val = DINT(DBLE(iwpara)/5.)
      if (MOD(iwpara,5).eq.0) then
          k = IDINT(val)
      else
          k = IDINT(val) + 1
      endif
      do i = 1,k
         imin = (i-1)*5 + 1
         imax = imin + 4
         if (imax.gt.iwpara) imax = iwpara
        write(1,'(A,1x,4(A," , "),A)') "# ",
     >  (wpara(j)(1:lnblnk(wpara(j))),j=imin,imax)
      enddo
       write(1,'(A)')"##############################################"   
c
c next it will use SBlim in VEGA 
       do i = 1, nfilt 
         mulim(i)=mulim(i) - abcor(i) 
       enddo
c
c GENERATING THE OUTPUT SIMULATED CATALOG
c
      write(6,*) 'simulated steps with mag:',nmagmax,' with z: ',nzmax
      nobj=0
      do i = 1,nmagmax             ! loop on magnitude in AB for iref
        mags= min + (i-1)*dmag
        magso = mags               ! mag AB iref 
        do j = 1, nzmax            ! loop on z 
           if (j.le.nztrans) then 
              zref=(j-1)*dz
              zbin = dz 
           else
              zref=(nztrans-1)*dz + dzsup*(j-nztrans)
              zbin = dzsup 
           endif   
           do k = 1, nlf           ! loop on LF type
              lf=lftyp(1,k)
              if (magso.lt.lftypm(1,k) .or.
     >            magso.gt.lftypm(2,k)) goto 23  
              if (zref.lt.lftypz(1,k) .or.
     >            zref.gt.lftypz(2,k)) goto 23  
c check the number of SEDs usable by LF type
              ntyp=0
              do p = lftyp(2,k), lftyp(3,k)
                if (zfor(linemod(1,p)).gt.0 .and.
     >              zfor(linemod(1,p)).gt.zref) then 
                     ntyp=ntyp+1
c                     if (lf==3)  write(*,*) p,linemod(1,p),linemod(2,p),
c     >  zfor(linemod(1,p)),zref,ntyp
                elseif (zfor(linemod(1,p)).le.0.) then
                     ntyp=ntyp+1
                endif   
              enddo  
              if (ntyp.eq.0) goto 23
c
              do p = lftyp(2,k), lftyp(3,k)          ! loop on model in given LF             
c                if (zfor(linemod(1,p)).lt.0 .or.
c     >              zfor(linemod(1,p)).le.zref)  goto 221
                if (zfor(linemod(1,p)).gt.0 .and.
     >              zfor(linemod(1,p)).le.zref)  goto 221
                 
                nexp = 0.
                nexti = 0
c check for various extinction E(B-V)
                do l = linemod(1,p), linemod(2,p) 
                   if (model(l).eq.p.and.zref.ge.z(l)
     >                  .and.zref.lt.z(l+1)) then
                      nexti = nexti + 1
                      lexti(nexti) = l
                   endif
                enddo
c take the mean  extinction for the density  
                if (nexti .eq. 1) then
                   lref=lexti(1)                 ! ident of ref model 
                else
c                   vali=IDNINT(DREAL(ran1(idummy))*DBLE(nexti))
                   if (MOD(nexti,2).eq.0) then 
                      vali = nexti/2
                   else
                      vali = (nexti-1)/2 + 1 
                   endif   
                   if (vali.ge.1 .and. vali.le.nexti) then 
                      lref=lexti(vali)           ! ident of ref model with exti
                   else
                      lref =lexti(1)
                   endif   
                endif   
                kcref =kcor(lfref,lref)          ! kcor lfref  and ref model
                kciref=kcor(iref,lref)           ! kcor iref  and ref model
                color =mag(lfref,lref)-mag(iref,lref) ! color AB=LFref-Iref
                ageref = age(lref)
                ebvref = exti(lref)
                zforref= zfor(lref)
c  expected density of object /deg^2/dz/dmag/type
c   LF given in LFRef AB  -> converting mag_iref_AB to lfref AB
c                mags=magso+color-kcref          ! mag AB iref -> lfref
cccc old version : Mabs(lfref)=mag(iref)+[mag(lfref)-mag(iref)]-kcor(lfref)-DM
c                            =mag(lfref)-kcor(lfref) -DM
c                 But this may be undefined due to the k-correction in lfref
cccc new version : Mabs(lfref)=mag(iref)-kcor(iref)+[mag(lfref)-mag(iref)]z=0 -DM
                mags=magso-kciref+colorz0(lfref,lref)  
c density 
                nexpd=surf*exp_dens(h0,om0,oml,zref,zbin,mags,
     >                 dmag,lf,mst0,phi0,alpha0,
     >                 lfevt,lfevz,lfphi,lfmst,lfalp,nlfev,
     >                 etal,etap,etas,
     >                 ecor,phicor,alpcor)
c                if (ntyp.gt.0) then
                   nexpd = nexpd/DBLE(ntyp)
c                else
c                   nexpd =0.
c                endif   
c  Number randomisation based on Poisson distribution with mean=nexp 
                if (nexpd.gt.1.e-6) then 
                  nexpo=SNGL(nexpd)
                  nexp= DBLE(poidev(nexpo,idummy))
                endif
c                write(6,*) magso,zref,nlf,k,lf,ntyp,p,lref,nexpd,nexp
                if (nexp.ge.1) then
                 do n = 1, IDNINT(nexp)                
                   rnum=p
                   typ=lfprof(lf)
c                  write(*,*) lf,typ
 10                dl= x(om0,oml,h0,zref)*(1+zref)
                   da=dl/(1+zref)**2
c  Uniform distribution of  AB iref mag within interval dmag
                   mags= magso+(ran1(idummy)-0.5)*dmag    ! mag AB iref  inside dmag  
c randomize Extinction for this model 
                   if (nexti.gt.1) then 
                     vali=IDNINT(DREAL(ran1(idummy))*DBLE(nexti))
                     if (vali.ge.1 .and. vali.le.nexti) then 
                       lref=lexti(vali)                ! ident of ref model with exti
                     else
                       lref =lexti(1)
                     endif   
                     kcref =kcor(lfref,lref)               ! kcor lfref  and ref model
                     kciref=kcor(iref,lref)                ! kcor iref  and ref model
                     color =mag(lfref,lref)-mag(iref,lref) ! color AB=LFref-Iref
                     ageref = age(lref)
                     ebvref = exti(lref)
                     zforref= zfor(lref)
                   endif   
c
c scale theoretical model to apparent magnitude 
                   scale=mags-mag(iref,lref)
c  computes the mag abs in Iref band 
                   mabs= mags-5*dlog10(dl)-25.-kciref ! Mabs(iref)_AB kecorrected
                   if (mabs.ge.0) then
                    goto 21
                   endif   
cccccccccccccccccccccccccccccccccc
c   extract all magnitudes in AB and put upper-limits 
                   do l = 1, nfilt
c  randomize AB magnitude within sigma_mag
                     magref(l) = mag(l,lref)+scale
                     magrefo(l)= magref(l)
                     mabsz(l)  = mabs + colorz0(l,lref)  ! Mabs(band)_AB kcorrected 
                     if (magref(l).le.uplim1(l)) then
                       temp=magref(l)
                       errmag(l) = 1.086/sig_noise(temp,l,Ae,Be,Ce,De)
                       magref(l) = error(l,magref,Ae,Be,Ce,De) 
                       if (magref(l).ge.uplim1(l)) then
                         temp=uplim(l)   
                         magref(l)=uplim(l) 
                         errmag(l) = -1                         
                       endif
                     else
                        magref(l)=uplim(l)
                        errmag(l) = -1  
                     endif  
cc   error and flux randomized 
                    if (magrefo(l).le.35) then 
                       temp=magrefo(l)
                    else
                       temp=35
                    endif  
                    bidf(l)=temp 
                    flran(l)=10**(-0.4*(temp+48.59))
                    errflran(l)= flran(l)/sig_noise(temp,l,Ae,Be,Ce,De)
                    flran(l) = error_flux(l,bidf,Ae,Be,Ce,De)
                   enddo   
cccccccccccccccccccccccccccccccccccccccccccccccccc
c                  
                   if (sb(1:1).eq.'y' .or. sb(1:1).eq.'Y') then
c  compute the re(pc) according to Mabs and type in Bj band 
                     mabs_bv= mabs-abcor(iref)+brcol(lref) ! Babs_veg from Iref
c  re0 in pc
                     re0=re_mabs(typ,zref,mabs_bv,h0,zeta)
c  re0as in arcsec
                     re0as=re0/(da*1.e6)*fac 
c   extract all magnitudes in AB with SB included
                     pm = 0 
                    do pmp = 1, 9
                       if (pmp.eq.1) then
                          mag_typ='ISO'
                          mag_fac=1.0
                          pm = pm + 1 
                       elseif (pmp.eq.2.and.compmag(2).eq.1) then
                          mag_typ='ISO'
                          mag_fac=rad_scale(1)
                          pm = pm + 1
                       elseif (pmp.eq.2.and.compmag(2).eq.0) then
                          goto 24 
                       elseif (pmp.eq.3.and.compmag(3).eq.1) then
                          mag_typ='AUTO'
                          mag_fac=rad_scale(2)
                          pm = pm + 1
                       elseif (pmp.eq.3.and.compmag(3).eq.0) then
                          goto 24 
                       elseif (pmp.eq.4 .and. compmag(4).eq.1) then
                          mag_typ='PETROSIAN'
                          mag_fac=eta_p
                          pm = pm + 1
                       elseif (pmp.eq.4 .and. compmag(4).eq.0) then
                          goto 24 
                       elseif (pmp.eq.5 .and. compmag(5).eq.1) then
                          mag_typ='APER'
                          mag_fac=rad_aper(1)
                          pm = pm + 1
                       elseif (pmp.eq.5 .and. compmag(5).eq.0) then
                          goto 24 
                       elseif (pmp.eq.6 .and. compmag(6).eq.1) then
                          mag_typ='APER'
                          mag_fac=rad_aper(2)
                          pm = pm + 1
                       elseif (pmp.eq.6 .and. compmag(6).eq.0) then
                          goto 24 
                       elseif (pmp.eq.7 .and. compmag(7).eq.1) then
                          mag_typ='APER'
                          mag_fac=rad_aper(3)
                          pm = pm + 1
                       elseif (pmp.eq.7 .and. compmag(7).eq.0) then
                          goto 24 
                       elseif (pmp.eq.8 .and. compmag(8).eq.1) then
                          mag_typ='APER'
                          mag_fac=rad_aper(4)
                          pm = pm + 1
                       elseif (pmp.eq.8 .and. compmag(8).eq.0) then
                          goto 24 
                       elseif (pmp.eq.9 .and. compmag(9).eq.1) then
                          mag_typ='APER'
                          mag_fac=rad_aper(5)
                          pm = pm + 1
                       elseif (pmp.eq.9 .and. compmag(9).eq.0) then
                          goto 24 
                       endif   
c
                     do l = 1, nfilt
c   --> get sigmat -> extract gt, gti from sigmat
                       sigmat=seeing(l)/(2.35*(re0/(da*1.e6)))/fac
                       if (sigmat.gt.sigtmax(1)) then
                          write(*,*) l,'Sigma_T >',sigtmax(1),':',sigmat
                          goto 21 
                       endif   
                       nsig=IDNINT(sigmat*10) + 1
c                 write(*,*) l,sigmat,sigt(1,nsig,1),sigt(1,nsig,ibmax)
                       do m = 1, ibmax
                          if (typ.eq.1) then
                            sigtr(m)=sigt(1,nsig,m)
                            lgtr(m)=lgt(1,nsig,m)
                            gtir(m)=gti(1,nsig,m)
                            betar(m)=beta(1,nsig,m)
                          else
                            sigtr(m)=sigt(2,nsig,m)
                            lgtr(m)=lgt(2,nsig,m)
                            gtir(m)=gti(2,nsig,m)
                            betar(m)=beta(2,nsig,m)
                          endif
                       enddo
c                       write(*,*) sigmat,nsig,sigtr(1),sigtr(ibmax)
                       mabs_bd=magrefo(l)-5*dlog10(dl)-25       ! Mabs(band)@z=0
     >                         -kcor(l,lref)  
                           ke=kcor(l,lref)    
                           sblim=mulim(l)+abcor(l)              ! SB_lim (AB) 
                           dmin=seeing(l)*fdmin(l)
                           call prof(typ,zref,mabs_bd,re0,da,ke,
     >                               sigtr,betar,lgtr,gtir,ibmax,
     >                               dmin,sblim,mag_fac,
     >                               mag_typ,h0,om0,oml,magp,
     >                               thetal,theta_mag) 
                           if (pmp.eq.1) theta(l) = thetal      ! diameter("): theta_iso
                           thetamag(pmp,l) = theta_mag          ! diameter("): theta_for mag
c upper-limits for ISO-PSEUDO-ISO and AUTO
                           if (pmp.le.4) then 
c  upper limits 
                              if (magp.gt.90) then
                                 magsb(pmp,l)=uplim(l)  
                                 errmagsb(pmp,l) = -1
                              else 
                                 magsb(pmp,l)= magrefo(l)+magp
c  adopting sextractor errors
                                 if (sexerr(1:1).eq.'Y'.or.
     >                                   sexerr(1:1).eq.'y') then 
                                   sblim = mulim(l)              ! SB_lim (Vega)
                                   call error_sim(thetamag,magsb,pmp,l,
     >            sblim,abcor,zp,gain,pixscale,fuge_fac,mag_ran,sn)
                                   magsb(pmp,l)=mag_ran
                                   errmagsb(pmp,l) = 1.086/sn
c adopting simple mag vs error behaviours 
                                 else
                                   magref1(l) = magrefo(l) + magp 
                                   temp=magref1(l)
                                   errmagsb(pmp,l) =
     >                              1.086/sig_noise(temp,l,Ae,Be,Ce,De)
                                   magsb(pmp,l) = error(l,magref1,
     >                                                Ae,Be,Ce,De) 
                                 endif
c last check for upper-limits 
                                 if (magsb(pmp,l).gt.uplim1(l)) then
                                    magsb(pmp,l)=uplim(l)  
                                    errmagsb(pmp,l) = -1
                                 endif
                              endif   
c upper-limits for APERTURES
                           elseif (pmp.gt.4) then
                              if (magp.gt.90) then
                                 magsb(pmp,l)=uplima(pmp-4,l)  
                                 errmagsb(pmp,l) = -1   
                              else
                                 magsb(pmp,l)= magrefo(l)+magp
c  adopting sextractor errors
                                 if (sexerr(1:1).eq.'Y'.or.
     >                                   sexerr(1:1).eq.'y') then 
c                                   
                                   sblim = mulim(l)              ! SB_lim (Vega)
                                   call error_sim(thetamag,magsb,pmp,l,
     >            sblim,abcor,zp,gain,pixscale,fuge_fac,mag_ran,sn)
                                   magsb(pmp,l)=mag_ran
                                   errmagsb(pmp,l) = 1.086/sn
c adopting simple mag vs error behaviours 
                                 else
                                   magref1(l) = magrefo(l) + magp 
                                   temp=magref1(l)
                                   errmagsb(pmp,l) =
     >                              1.086/sig_noise(temp,l,Ae,Be,Ce,De)
                                   magsb(pmp,l) = error(l,magref1,
     >                              Ae,Be,Ce,De) 
                                 endif
c last check for upper-limits 
                               if (magsb(pmp,l).gt.uplima1(pmp-4,l))then
                                    magsb(pmp,l)=uplima(pmp-4,l)  
                                    errmagsb(pmp,l) = -1
                               endif
                              endif                                    
                           endif   
                     enddo   ! l : nb of filters   
 24                  continue
                    enddo    ! pm: magnitude type
                  endif      ! sb : surface brightness

c              write(1,'(2(2(f8.3,1x),i5,1x),3(f8.3,1x),
c     >                 500(E12.5,1x))') 
c     >  mags,zref,lf,kciref,mabs,rnum,ecor,phicor,alpcor,
c     > (magrefo(l),l=1,nfilt),(mabsz(l),l=1,nfilt),
c     >  (magref(l),errmag(l),l=1,nfilt),
c     >  re0,mabs_bv,(thetamag(1,l),l=1,nfilt),theta_iso(l) 
c     > ((magsb(m,l),errmagsb(m,l),l=1,nfilt),m=1,nbmag)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c writing in output file 
c  test if object writes in output 
                  pm=0
                  do m=1,nfilt
                     if (errmag(m).le.0) pm=pm+1
                  enddo   
                  if (pm.eq.nfilt) then
c                     write(*,*) 'no object' 
                     goto 21 
                  endif
                  nobj=nobj + 1
                  call WRITE_SIM(wpara,iwpara,sb,nobj,
     >  mags,zref,lf,kciref,ageref,ebvref,zforref,
     >  mabs,rnum,ecor,phicor,alpcor,
     >  nfilt,magrefo,mabsz,re0,re0as,mabs_bv,thetamag,
     >  magref,errmag,flran,errflran,magsb,errmagsb,dl,
     >  str_out,iwout)
c
                  write(1,'(100(A,2x))') 
     >            (str_out(m)(1:lnblnk(str_out(m))),m=1,iwout)
cccccccccccccccccccccccccccccccccccccccccc
c screen info 
                  write(6,1010) nobj,mags,zref,lf,nexp,rnum,char(13)
                  call flush(6)
cccccccccccccccccccccccccccccccccccccccccc
 21               continue
               enddo    ! close loop (n) on Nb objects to be simulated 
             endif      ! only if n>=1
 221         continue
            enddo       ! close loop (p) on type in LF 
 23         continue
           enddo        ! close loop (k) on LF type
        enddo           ! close loop (j) on z 
      enddo             ! close loop (i) on magnitude 
      close(1)
      write(*,*) '                              ' 
      write(*,*) ' N object = ', nobj
      write(*,*) '  That s all folks !!!  '
c
 1010 format("object -> ",i8,1x,2(f8.3,1x),I6,1x,f8.3,1x,I6,a1,$)
 1011 format("object -> ",i8,1x,2(f8.3,1x),I6,1x,f8.3,1x,I6)
c
      stop
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine   ERROR_SIM(thetamag,magsb,pm,filt,sblim,abcor,zp,
     >                       gain,pixscale,fuge_fac,mag_ran,sn)
c  randomizing MAG_AB  according to errors including 
c     area sigma_bg, gain and flux
      implicit none 
      integer*4 nbf,pm,filt
      INCLUDE   'dim_filt.decl'
      real*8  mag_ran,sn,err,area,sigbg,flux,fluxab,dflux,fluxn
      real*8  zp(nbf),gain(nbf),pixscale(nbf)
      real*8  fuge_fac(nbf)
      real*8  abcor(nbf),thetamag(10,nbf),magsb(10,nbf)
      real*8  sblim,pi,errmin
      real     gasdev
      external gasdev
c
      errmin=0.03
      pi=3.141592654      
      area=pi*(thetamag(pm,filt)/2/pixscale(filt))**2       !diam(")->area(pix)
      sigbg=10**(-0.4*(sblim-zp(filt)-5*dlog10(pixscale(filt)))) !sig_bg(adu)
      sigbg=sigbg*fuge_fac(filt)                            ! sig_bg*fuge_factor
      flux=10**(-0.4*(magsb(pm,filt)-abcor(filt)-zp(filt))) !flux(adu)
      err=1.0857/flux*dsqrt(area*sigbg**2+flux/gain(filt))  !error
      err = dsqrt(err**2+errmin**2)                         !adds errors in quad.
      sn = 1.0857/err                                       !signal to noise
c   randomize mag according to sigma 
      fluxab=10**(-0.4*(magsb(pm,filt)+48.59))              !flux (AB)
      dflux= fluxab*err/1.0857                              !err_flux
 13   fluxn = fluxab  + 1.0*gasdev(SNGL(dflux))             !flux+ran (AB) 
      if (fluxn.le.0) goto 13
      mag_ran=-2.5*dlog10(fluxn) - 48.59                    !Mag+ran (AB)
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine   MAGLIM_SIM(theta,thres,filt,sblim,abcor,zp,
     >                       gain,pixscale,fuge_fac,mag_lim)
c  computing MAG_LIM in AB  according to area sigma_bg, gain and flux
      implicit none 
      integer*4 nbf,filt,i
      INCLUDE   'dim_filt.decl'
      real*8  sn,err,area,sigbg,flux
      real*8  zp(nbf),gain(nbf),pixscale(nbf)
      real*8  fuge_fac(nbf),thres,theta
      real*8  abcor(nbf),mag,mag_lim
      real*8  sblim(nbf),pi,errmin
c
      errmin=0.03
c
      pi=3.141592654      
c
      area=pi*(theta/2/pixscale(filt))**2                      !diam(")->area(pix)
      sigbg=10**(-0.4*(sblim(filt)-zp(filt)-5*dlog10(pixscale(filt)))) !sig_bg(adu)
      sigbg=sigbg*fuge_fac(filt)                               !sig_bg*fuge_factor
      
      do i =1, 1000
         mag = 18 + i*0.02                                     ! mag AB
         flux=10**(-0.4*(mag-abcor(filt)-zp(filt)))            !flux(adu)
         err=1.0857/flux*dsqrt(area*sigbg**2+flux/gain(filt))  !error
         err = dsqrt(err**2+errmin**2)                         !adds errors in quad.
         sn = 1.0857/err                                       !signal to noise
         if (sn.le.thres) then
           mag_lim=mag                                         !Mag+ran (AB)
           goto 10
         endif  
      enddo
c
 10   return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL*8 FUNCTION  ERROR(n,mag,A,B,C,D)
c   return the magnitude + random error in AB 
c   inside gaussian distribution for the noise   
      implicit none 
      integer*4 n,nbf
      INCLUDE   'dim_filt.decl'
      real*8     mag(nbf)
      real*8     fab(nbf),dfab,magv
      real*8     flux
      real*8     sig_noise
      real*4       gasdev,bid,bid1
      real*8     A(nbf),B(nbf),C(nbf),D(nbf)
      external   sig_noise,gasdev

c  adding random errors to mag AB   
      fab(n)=10**(-0.4*(mag(n)+48.59))
      magv=mag(n)
      dfab=fab(n)/sig_noise(magv,n,A,B,C,D)
      bid1 = SNGL(dfab)
 13    bid=gasdev(bid1)
      flux = fab(n)  + bid
      if (flux.le.0) goto 13
c      write(*,*) 'sig:',bid
      error=-2.5*dlog10(flux) - 48.59
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL*8 FUNCTION  ERROR_FLUX(n,mag,A,B,C,D)
c   return the magnitude + random error in AB 
c   inside gaussian distribution for the noise   
      implicit none 
      integer*4 n,nbf
      INCLUDE   'dim_filt.decl'
      real*8     mag(nbf)
      real*8     fab(nbf),dfab,magv
c      real*8     flux
      real*8     sig_noise
      real*4       gasdev,bid,bid1
      real*8     A(nbf),B(nbf),C(nbf),D(nbf)
      external   sig_noise,gasdev

c  adding random errors to flux AB   
      fab(n)=10**(-0.4*(mag(n)+48.59))
      magv=mag(n)
      if ((sig_noise(magv,n,A,B,C,D)).le.0.2) then
        dfab=5*fab(n)
c        dfab=fab(n)/sig_noise(magv,n,A,B,C,D)
        fab(n) = 0.
      else
        dfab=fab(n)/sig_noise(magv,n,A,B,C,D)
      endif
      bid1 = SNGL(dfab)
 13   bid  = gasdev(bid1)
      error_flux = fab(n)  + bid
c      if (abs(dfab/fab(n)).le.10) goto 13
c      write(*,*) 'sig:',bid
c      error=-2.5*dlog10(flux) - 48.59
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL*8 FUNCTION  MAGLIM(band,thres,A,B,C,D)
c  computing upper limit according to the sigma limit 
      implicit none
      integer*4 band,i,nbf
      INCLUDE   'dim_filt.decl'
      real*8    sig_noise,mag,thres
      real*8    A(nbf),B(nbf),C(nbf),D(nbf)
      external sig_noise
      maglim=0
      mag=10
      do i = 1, 20000
          mag=mag+0.01
          if (sig_noise(mag,band,A,B,C,D).le.thres) then
             maglim=mag
             return
          endif
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL*8 FUNCTION SIG_NOISE(magn,band,A,B,C,D)
c  compute err as a function of magnitude assuming err=schechter function
c  errors are given in AB magnitude 
      implicit none
      integer*4 band,nbf  
      INCLUDE   'dim_filt.decl'
      real*8 magn,A(nbf),B(nbf),C(nbf),D(nbf)
      real*8 errmin,temp
c
c computing mag_err assuming mag_err=schechter function (A,B,C,D)
c   and S/N = 1.0857/mag_err
      errmin=0.03
      if (magn.le.A(band)) then 
        temp = 10**(0.4*(C(band)+1)*(magn-A(band)))
        sig_noise=B(band)*10**(0.4*(C(band)+1)*(magn-A(band)))
      else
        temp = 10**(D(band)*(magn-A(band)))
        if(10**(D(band)*(magn-A(band))).le. 40.) then
         sig_noise=B(band)/2.71828*dexp(10**(D(band)*(magn-A(band))))
        else
         sig_noise=B(band)/2.71828*dexp((40.d0))
        endif
c        sig_noise=B(band)/2.71828*dexp(10**(D(band)*(magn-A(band))))
      endif  
      sig_noise=dsqrt(sig_noise**2+errmin**2)
c      if (sig_noise.le.errmin) sig_noise=errmin
      sig_noise=1.0857/sig_noise 
c
      return
      end         
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL*8 FUNCTION exp_dens(h0,om0,oml,z,dz,mag,
     >       dmag,lf,mst0,phi0,alpha0,
     >       lfevt,lfevz,lfphi,lfmst,lfalp,nlfev,
     >       etal,etap,etas,
     >       ecor,phicor,alpcor)
c
c   estimate of the number of gal/dmag/dz/type/deg^2 based on LF 
      implicit none 
      integer*4 i,lf,lfevt(100),nlfev
      real*8    alpha,mstar,mabs,dmabs,mag,dmag
      real*8    phi,phistar
      real*8    c,h0,om0,oml,x,dl,da,z,dz,dvz
      real*8    ecor,phicor,alpcor,etal,etap,etas
      real*8    lfevz(100),lfphi(100),lfmst(100),lfalp(100)
      real*8    mst0(100),phi0(100),alpha0(100)
      real*8    mstar0,mstar1,phistar0,phistar1,alp0,alp1,z0,z1
      external  x
c
      c=300000.0
      exp_dens=0.
      alpha   =0.
      phistar =0.
      mstar   =0.
      ecor    =0.
      phicor  =1.
      alpcor  =1.
        z0      = 0.0
        mstar0  =  0
        phistar0=  0
        alp0  =  0
        z1      =  0
        mstar1  =  0
        phistar1=  0
        alp1  =  0
c
      if (nlfev .eq. 0 ) then
c L*(z)=L*(0).(1+z)^-etal  -> M*(z) = M*(0) + 2.5.etal*log(1+z) 
         ecor    = 2.5*etal*dlog10(1+z)
c         mstar   = mst0(lf)   + 5*dlog10(h0/100) + ecor 
         mstar   = mst0(lf)   + ecor 
c Phi*(z)=Phi*(0)*(1+z)^etap
         phicor  = (1+z)**etap
c         phistar = phi0(lf) * (h0/100)**3 * phicor
         phistar = phi0(lf) * phicor
c alpha(z) = alpha(0) + etas*log(1+z)
         alpcor  = etas*dlog10(1+z)
         alpha   = alpha0(lf) + alpcor
      else   
        z0      = 0.0
        mstar0  =  mst0(lf) 
        phistar0=  phi0(lf) 
        alp0  =  alpha0(lf)
        z1      = 99.0
        mstar1  =  mst0(lf) 
        phistar1=  phi0(lf) 
        alp1  =  alpha0(lf)
       do i = 1, nlfev
         if (lfevt(i).eq.lf .and.
     >          z.ge.lfevz(i) ) then
            z0      =  lfevz(i)
            mstar0  =  lfmst(i) 
            phistar0=  lfphi(i) 
            alp0    =  lfalp(i)
         endif
         if (lfevt(i).eq.lf .and.
     >          z.lt.lfevz(i) .and. lfevz(i).gt.z0) then
            z1      =  lfevz(i)
            mstar1  =  lfmst(i) 
            phistar1=  lfphi(i) 
            alp1    =  lfalp(i)
            goto 37
         endif     
       enddo
 37    continue
       mstar  = mstar0   + (mstar1-mstar0)/(z1-z0)*(z-z0)
       phistar= phistar0 + (phistar1-phistar0)/(z1-z0)*(z-z0)
       alpha  = alp0     + (alp1-alp0)/(z1-z0)*(z-z0)
c       mstar  = mstar    + 5*dlog10(h0/100)  
c       phistar= phistar  * (h0/100)**3 
       alpcor = alpha
       phicor = phistar
       ecor   = mstar
      endif
c
c   computing distance Luminosity and angular 
      da= x(om0,oml,h0,z)/(1+z)
      dl= x(om0,oml,h0,z)*(1+z)
c
c   computing the comoving volume element : dV/dz   
      dvz=1./dsqrt(1+om0*z-oml+oml/(1+z)**2)
      dvz = (c/h0)*da*da*(1+z)*dvz         
c
c  computing Mabs
      mabs  = mag -5*dlog10(dl) -25.
c    
      dmabs = mstar-mabs
c
c  computing  phi
      phi=dexp(-10**(0.4*dmabs))*10**(0.4*dmabs*(alpha+1))
      phi= 0.4*dlog(10.d0)*phi*phistar
c  expected number /dmag/dz in area square deg
      exp_dens=phi*dvz*dz*dmag/3282.8064 
c      write(*,'(4(f6.2,1x),2(E12.6,1x),f6.3,1x,E12.6,1x,f8.1)')
c     >  z,dz,mag,dmag,dvz,dl,kcor,phi,exp_dens
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccc
      REAL*8 function re_mabs(t,z,m,h0,zeta)
c  relation between re-Mabs(B band)
c  -Mb = p log(re(kpc)) +q + (p-5)*log(Ho/50)
c
      implicit none
      integer*4 t
      real*8    m,p,q,h0,z,zeta,sigmuo,ms
      real     gasdev
      external gasdev 
c
      p=0.
      q=0.
      sigmuo=0.35
c

      if (t.eq.4) then        ! For spheroids based on Binggeli et al. (1984) 
         ms=-20.5+5*dlog10(h0/50)
         if (m.ge.ms) then ! Ms<-20.5, Bertin et Arnouts (AASS 117, 1996) (A3)
            p=10
            q=15.5
         else                 ! Ms<-20.5, Bertin et Arnouts (AASS 117, 1996) (A3)
            p=3.33
            q=18.833
         endif   
        re_mabs=1000*10**(-1./p*(q+m+(p-5)*alog10(SNGL(h0)/50.))) 
         ms=-23.+5*dlog10(h0/50)
        if (m.lt.ms) then    ! maximal value of re if Mb<-23 (->Mb=-23) 
          re_mabs=1000*10**(-1./p*(q-23+(p-5)*alog10(SNGL(h0)/50)))
        endif  
      else                    ! For disk based on Freeman law (1970) muo=21.65
c                               muo=Mabs+5log(re[pc]) + 22.44
c                               log(re(pc)) =-0.2(basb+0.79) 
c                               log(re(kpc))=-0.2*(Babs+15.79)
c add contribution from dwarves if Mb>-17 with 
c  muo = 21.65 + 0.7*(Mb+15)
         ms=-17+5*dlog10(h0/50)
         if (m.le.ms) then  ! Freeman non dwarves 
            p=5
            q=15.79
         elseif (m.gt.ms) then  ! Freeman dwarves 
            p=16.667
            q=12.97     ! if Mabs>-17
c            q=17.63     ! if Mabs>-15
         endif   
         re_mabs=1000*10**(-1./p*(q+m+(p-5)*alog10(SNGL(h0)/50)))
         re_mabs=re_mabs*(1+z)**zeta      
      endif
 10   re_mabs=re_mabs*(1 + gasdev(SNGL(sigmuo))*0.4605)
c      write(*,*) re_mabs
      if (re_mabs.le.5) goto 10 
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION gasdev(sig)
      INTEGER idummy
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1,sig
      SAVE iset,gset
      DATA iset/0/
      common /dum/        idummy
      if (iset.eq.0) then
1       v1=2.*ran1(idummy)-1.
        v2=2.*ran1(idummy)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      gasdev=gasdev*sig
c
      return
      END
c
ccccccccccccccccccccccccccccccccccccc
        real*8 function x(omo,oml,ho,z)
ccccccccccccccccccccccccccccccccccccc
        implicit none
        integer*4 i 
        real*8 omo,ho,z,c,oml,omt,ao,sum,dz,zi,Ez
c       
        x=0
        c = 300000.
        omt = oml + omo
        if (z.le.1.e-5) then 
           x=1.e-5
        else   
           if (omt.lt.1.and.oml.eq.0) then
c          ao = c/(ho*sqrt(ABS(1-omt)))
c  in fact we use x = ao * x(z) with x(z) from eq 8 of 
c    Moscardini et al.  So we don't need to compute ao   
             ao =1
             x = (omo*z-(omo-2)*(1-dsqrt(1+omo*z)))/(omo*omo*(1+z))
             x= 2*c/(ao*ho)* x
           elseif (DABS(omt-1).le.1.e-5.and.oml.eq.0) then
             ao = 1
             x = (omo*z-(omo-2)*(1-sqrt(1+omo*z)))/(omo*omo*(1+z))
             x= 2*c/(ao*ho)* x
           elseif (omo.le.1.and.oml.ne.0) then
             ao  = 1
             sum = 0 
             dz  = z/50. 
             do i = 1,51
               zi = (i-0.5)*dz
               if (zi.gt.z) goto 1                     
               Ez=dsqrt(omo*(1+zi)**3+(1-omo-oml)*(1+zi)**2+oml)
               sum = sum + dz/Ez 
             enddo        
 1           x =  c/(ho*ao) * sum 
           endif
        endif   
c         
        return
        end
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION poidev(xm,idum)
      INTEGER idum
      REAL poidev,xm,PI
      PARAMETER (PI=3.141592654)
CU    USES gammln,ran1
      REAL alxm,em,g,oldm,sq,t,y,gammln,ran1
      SAVE alxm,g,oldm,sq
      DATA oldm /-1./
      if (xm.lt.12.)then
        if (xm.ne.oldm) then
          oldm=xm
          g=exp(-xm)
        endif
        em=-1
        t=1.
2       em=em+1.
        t=t*ran1(idum)
        if (t.gt.g) goto 2
      else
        if (xm.ne.oldm) then
          oldm=xm
          sq=sqrt(2.*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.)
        endif
1       y=tan(PI*ran1(idum))
        em=sq*y+xm
        if (em.lt.0.) goto 1
        em=int(em)
        t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
        if (ran1(idum).gt.t) goto 1
      endif
      poidev=em
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
c
ccccccccccccccccccccccccccccccccccccccc
      FUNCTION factrl(n)
      INTEGER n
      REAL factrl
CU    USES gammln
      INTEGER j,ntop
      REAL a(33),gammln
      SAVE ntop,a
      DATA ntop,a(1)/0,1./
      factrl=1
      if (n.lt.0) then
        pause 'negative factorial in factrl'
      else if (n.le.ntop) then
        factrl=a(n+1)
      else if (n.le.32) then
        do 11 j=ntop+1,n
          a(j+1)=j*a(j)
11      continue
        ntop=n
        factrl=a(n+1)
      else
        factrl=exp(gammln(n+1.))
      endif
      return
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
c
ccccccccccccccccccccccccccccccccccccc
      subroutine prof(t,z,mabs,re,dao,kcor,
     >  sig,beta,lgt,gti,imax,
     >  dmin,sblim,mag_fac,
     >  mag_typ,ho,omo,oml,magc,thetalim,theta_mag)
      implicit none 
      integer*4  t,max,i,j,imax,m
      parameter (max=310)
      real*8     mabs,kcor,z,sblim,dmin,mag_fac,magc
      real*8     sig(max),beta(max),lgt(max)
      real*8     gti(max),gtinf,betamag
      real*8     x,dao,re,ho,omo,oml
      real*8     pi,a,lre,rk_1,gt_int,theta_int,beta_1,eta
      real*8     muz,thetalim,theta_mag,gtr_int,dbeta,min_fac
      real*4     factrl
      character  mag_typ*512

      external   x,factrl
c
      pi=3.141592654
      a=0
      beta_1=0
c  compute Gt(infinity)
      if (t.eq.1) then
         m=2
         a=1.68
      elseif (t.eq.4) then
         m=8
         a=7.67
      endif   
      gtinf=factrl(m)*pi/(a**m)
       thetalim=-99
       magc=0
c       write(6,*) t,z,mabs,re,kcor
c  test if objects always detected 
       lre=dlog10(re)
c
         muz = mabs +kcor +5*dlog10(re) +10*dlog10(1+z) 
     >          + 21.5721 +2.5*dlog10(gtinf) -2.5*lgt(imax)
         if (muz.le.sblim) then
            magc=0
            thetalim=999
            goto 10
         endif
c
       do i = 1, imax
          muz = mabs +kcor +5*dlog10(re) +10*dlog10(1+z) 
     >          +21.5721 +2.5*dlog10(gtinf) -2.5*lgt(i)
          if (i.eq.1 .and. muz.gt.sblim) then 
             thetalim=0
             theta_mag=0
             magc = 99
             goto 10
          endif   
          if (muz.ge.sblim) then
             thetalim=2*(re/(dao*1.e6))*beta(i)*180*3600/pi ! diam(") 
             theta_mag=thetalim                             ! diam(")
             if (thetalim.lt.dmin) then                     ! seeing diameter
                  magc=99
                  goto 10 
             endif     
          endif   
c   Pseudo-ISOPHOTAL MAGNITUDES 
          if (mag_typ(1:3).eq.'ISO' .or. mag_typ(1:3).eq.'iso') then  
c
             if (muz.ge.sblim) then
               thetalim=2*(re/(dao*1.e6))*beta(i)*180*3600/pi ! diam(")
               theta_mag = thetalim*mag_fac                   ! diam(")
               betamag=mag_fac*beta(i)                        ! radius(beta)
               if (mag_fac.eq.1) then
                  magc=-2.5*dlog10(gti(i))
                  goto 10 
               elseif (mag_fac.lt.1) then
                  do j = 1,i
                    if (beta(j).ge.betamag) then
                       magc=-2.5*dlog10(gti(j))
                       goto 10 
                    elseif (j.eq.imax.and.beta(j).lt.betamag) then
                       magc=-2.5*dlog10(gti(imax))
                       goto 10 
                    endif
                  enddo
               elseif (mag_fac.gt.1) then
                  do j = i,imax
                     if (beta(j).ge.betamag) then
                        magc=-2.5*dlog10(gti(j))
                        goto 10 
                     elseif (j.eq.imax.and.beta(j).lt.betamag) then
                        magc=-2.5*dlog10(gti(imax))
                        goto 10
                     endif
                  enddo
               endif
             endif   
c   APERTURE MAGNITUDES --> mag_fac= diameter aperture in arcsec
         elseif (mag_typ(1:4).eq.'APER'
     >           .or. mag_typ(1:4).eq.'aper') then
c
             betamag=mag_fac/2                ! radius(")
             theta_mag= betamag*2             ! diam(")
             betamag=mag_fac/2/(re/(dao*1.e6))/(180*3600/pi) !radius(beta)
             if (beta(i).ge.betamag) then
                magc=-2.5*dlog10(gti(i))
                goto 10 
             elseif (i.eq.imax.and.beta(i).lt.betamag) then
                magc=-2.5*dlog10(gti(imax))
                goto 10
             endif
c    AUTOMATIC  MAGNITUDES                    
         elseif (mag_typ(1:4).eq.'AUTO'
     >                  .or. mag_typ(1:4).eq.'auto') then
c  
            min_fac=1.5 
            if (muz.ge.sblim) then
c get theta_lim= 2.5xtheta_iso for integration of first moment rk 
               thetalim=2.5*(re/(dao*1.e6))*beta(i)*180*3600/pi  !radius(")
c compute the first moment rk   
               gt_int=0
               gtr_int=0
               dbeta=beta(2)-beta(1)
               do j = 1, imax
c  integrating g~(theta) up to 2.5*theta_iso
                  gt_int= gt_int + 10**(lgt(j))*dbeta
                  gtr_int= gtr_int + 10**(lgt(j))*beta(j)*dbeta
c  get integral G~(theta)= int(0 to 2*theta_iso) theta * g~(theta)  / 2pi 
c  and first moment radius 
                  theta_int = (re/(dao*1.e6))*beta(j)*180*3600/pi !radius(") 
                  if (theta_int.ge.thetalim) then
c                     rk_1 = gti(j)/(2*pi)*gtinf / gt_int 
                     rk_1 = gtr_int / gt_int 
                     beta_1 = rk_1                           ! radius(beta)
                     rk_1 = rk_1*(re/(dao*1.e6))*180*3600/pi ! radius(")
                     theta_mag = rk_1*mag_fac                ! radius(") 
                     goto 25 
                  endif   
               enddo    
 25            betamag=mag_fac*beta_1      ! radius(beta) 
               thetalim=thetalim/2.5       ! radius_iso(") 
c if theta < min_facxtheta_iso -> used min_facxtheta_iso
               if (theta_mag .le. (min_fac*thetalim)) then
                  betamag=min_fac*thetalim/(re/(dao*1.e6)*180*3600/pi)
                  theta_mag=2*min_fac*thetalim    ! diam(")
                  thetalim=2*thetalim             ! diam(")
               else
                  theta_mag=2*theta_mag           ! diam(")
                  thetalim=2*thetalim             ! diam(")                  
               endif   
               do j = 1, imax
                 if (beta(j).ge.betamag) then
                   magc=-2.5*dlog10(gti(j))
                   goto 10 
                 elseif (j.eq.imax.and.beta(j).lt.betamag) then
                   magc=-2.5*dlog10(gti(imax))
                   goto 10
                 endif      
               enddo
            endif   
c Petrosian magnitudes
         elseif (mag_typ(1:5).eq.'PETRO'
     >                  .or. mag_typ(1:5).eq.'petro') then 
            min_fac=2
            if (i.gt.1) then 
               eta=10**(lgt(i))/(gti(i)*gtinf)*pi*beta(i)**2
            else   
               eta=999 
            endif   
            if (eta .le. mag_fac) then 
              theta_mag=2*min_fac*beta(i)*(re/(dao*1.e6))*180*3600/pi
              betamag=min_fac*beta(i)
c              write(6,'(3(f6.3,1x))') eta,betamag,theta_mag
               do j = 1, imax
                 if (beta(j).ge.betamag) then
                   magc=-2.5*dlog10(gti(j))
                   goto 10 
                 elseif (j.eq.imax.and.beta(j).lt.betamag) then
                   magc=-2.5*dlog10(gti(imax))
                   goto 10
                 endif      
               enddo   
            endif   
         endif
       enddo
 10   return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   READ GALAXY LIBRARY 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_libsim(lib,size,zstep,zmax,dz,
     > hub,omeg,lbda,
     > model,zlib,maglib,klib,zfor,exti,age,
     > namefl,magtyp,imagm,recmax)
c
      implicit none 
      integer*4     lnblnk,i,k,size,test,nmodel,nzf
      integer*4     imagm,vali,reclibmax
      character*512 zpdir,zpwork,file,valc
      integer*4     nbf
      INCLUDE 'dim_filt.decl'      
      integer*4     recmax,modext(2),nebv
      real*8        zstep,zmax,dz,hub,omeg,lbda,ebv(100),zf(100)
      character*512 lib,extlaw
      character*512 magtyp,namefl(nbf)
      character*512 paravc(100),param
      integer*4     model(size)
      real*8        zfor(size),exti(size)
      real*8        zlib(size),maglib(nbf,size)
      real*8        klib(nbf,size),age(size)
c
c Environemental  Variable 
      call getenv('LEPHAREDIR',zpdir)
      if (lnblnk(zpdir) .eq. 0) then
        write(6,*) 'WARNING :  variable LEPHAREDIR not defined'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(6,*) 'WARNING :  variable LEPHAREWORK not defined'
        stop
      endif
c
c reading doc files 
        file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >            // lib(1:lnblnk(lib)) // '.doc'
c
        param='LIB_TYPE'
        call read_para2(param,file,paravc,test)
        valc=paravc(1)       
        if (valc(1:3).eq.'GAL') then
          param='NUMBER_SED'
          call read_para(param,file,paravc,test)
          read(paravc(1),'(i8)') nmodel  
          param='NUMBER_ROWS'
          call read_para(param,file,paravc,test)
          read(paravc(1),'(i8)') recmax
          param='RECORD_LENGTH'
          call read_para(param,file,paravc,test)
          read(paravc(1),'(i8)') reclibmax 
          param='Z_STEP'
          call read_para2(param,file,paravc,test)
          read(paravc(1),'(f8.3)') zstep
          read(paravc(2),'(f8.3)') zmax
          read(paravc(3),'(f8.3)') dz
          param='COSMOLOGY'
          call read_para2(param,file,paravc,test)
          read(paravc(1),'(f8.3)') hub
          read(paravc(2),'(f8.3)') omeg
          read(paravc(3),'(f8.3)') lbda
          param='Z_FORM'
          call read_para2(param,file,paravc,test)
          nzf=test
          if (nzf.ge.1) then
            do i = 1,nzf
               read(paravc(i),'(f8.3)') zf(i)
            enddo
          else
             nzf = 0 
             do i = 1, 100
                zf(i) = 0.
             enddo   
          endif
          param='EXTINC_LAW'
          call read_para2(param,file,paravc,test)
          extlaw = paravc(1) 
          param='EB_V'
          call read_para2(param,file,paravc,test)
          nebv=test
          do i = 1,nebv
             read(paravc(i),'(f8.3)') ebv(i)
          enddo
          param='MOD_EXTINC'
          call read_para2(param,file,paravc,test)
          read(paravc(1),'(i8)') modext(1)
          read(paravc(2),'(i8)') modext(2)
           param='MAG_TYPE'
           call read_para2(param,file,paravc,test)
           magtyp=paravc(1)
           param='FILTERS'
           call read_para2(param,file,paravc,test)
           imagm=test
           do i = 1,test
              namefl(i) = paravc(i)
           enddo
        else
           write(*,*) ' NOT A GALAXY s LIBRARY --> STOP '
           stop
        endif 
        close(1)
        if (recmax.gt.size) then
          write(6,*) ' dim lib :',recmax,' default size :',size
          write(6,*) '-> The dimension of library exceeds the'
          write(6,*) ' default s one  in dim_lib.decl  '
          write(6,*) ' dim_lib.decl must be increased  '
          STOP
        endif   
c reading the library 
c
        file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     > // lib(1:lnblnk(lib)) // '.bin'
        open(2,file=file,status='unknown',access='direct',
     >               recl=reclibmax)
        do i = 1, recmax 
          if (nzf.ge.1) then 
            read(2,rec=i) model(i),exti(i),
     >                  zlib(i),age(i),zfor(i),vali,
     >                  (maglib(k,i),k=1,vali),
     >                  (klib(k,i),k=1,vali)
c            write(*,*) 'read zform ...'
          else
            read(2,rec=i) model(i),exti(i),
     >                  zlib(i),age(i),vali,
     >                  (maglib(k,i),k=1,vali),
     >                  (klib(k,i),k=1,vali)
            zfor(i) = 0. 
c            write(*,*) 'read no zform ...'
          endif   
        enddo
        close(2) 
        write(*,*) ' library length:',recmax
c
      RETURN
      END
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




