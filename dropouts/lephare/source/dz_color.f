c
c  last modif 11/09/07
c     Creation of Dz for a given SED and filter set and magnitudes
c     The program needs as input : 
c     -1- A magnitude library through all filters  
c     -2- An input SED type 
c     -3- The behaviour of errors versus magnitude 
c
      PROGRAM DZ_COLOR
      implicit none
c
      integer*4 nbf,inlib,sbmax,nlf
      INCLUDE   'dim_filt.decl'
      INCLUDE   'dim_lib.decl'
      integer   idummy
      parameter (sbmax=310)
      integer*4 i,j,j2,k
      integer*4 kmax,imagm                 
      integer*4 nfilt,lnblnk
      integer*4 model(inlib)
      real*8    mag(nbf,inlib),z(inlib),kcor(nbf,inlib)
      real*8    zfor(inlib),exti(inlib),age(inlib)
      real*8    Ae(nbf),Be(nbf),Ce(nbf),De(nbf)
      real*8    flmoy(nbf),flwidth(nbf)
      real*8    maglim
      real*8    val,error,sig_noise
      real*8    x,c,h0,om0,oml,dz,dzsup,zmod
      real*8    exp_dens
      real*8    hzsteidel,hz3lf
      real*8    fcorr(nbf)
c SB effect parameters
      real*8    re_mabs,fac,pi
c
c ERROR SIM parameters 
      real*8    zp(nbf),abcor(nbf)
c
c other parameters 
      real       poidev
      real*4     ran1
      real*8     schechter
c
      character*512  inf,simf,errf
      character*512  type
c
      real*8        paravr(100)  
      integer*4     paravi(100)
      character*512 paravc(100),namefl(nbf),nameff(nbf)
      character*1024 str

c
      character*512 zpdir,zpwork,config,param,filters
c
      integer*4  test,jmax
ccccccccccccccccc
      integer*4 sed_lev(2),filt_lev
      real*8    mag_lev,zmax_lev
      real*8    mag_temp,err_temp,err_temp2
      real*8    scale,dztot,dzcol,ecol,dzcol1,dzcol2
c
      common /dum/        idummy
      common /errorfile/  errf
      EXTERNAL schechter,sig_noise,maglim
      external x,exp_dens,error,poidev,hzsteidel,hz3lf,re_mabs,ran1
cccccccccccccccccccccccccccccccccc
c constant parameters
      c=300000.
      pi=3.141592654
      fac=180*3600/pi
c
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
      write(*,'(A)') config(1:lnblnk(config)) 

      write(6,'(2(A,1x))') '   >>> reading keywords from ',
     > config(1:lnblnk(config))
c
c
cccccccccccccccccccccccccccccccccccccc
c  Reading Filter FILE 
      param='-FILTER_FILE'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  filters=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
      write(6,'(A)') '   >>> filter characteristics ...'
c  Extract names of filters 
        str = zpwork(1:lnblnk(zpwork)) //'/filt/'
     >            // filters(1:lnblnk(filters)) 
        open(1,file=str,status='old')
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
      call zeropoint(filters,zp,abcor,flmoy,flwidth,fcorr,nfilt)
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
      param='-ERROR_FILE'
      call getc_option(param,config,100,paravc,test)
      if (test.eq.1 ) then
           errf=paravc(1)(1:lnblnk(paravc(1)))
           call get_path(errf)
      else
         call err_option(param,1)         
      endif     
c
      param='-SED_TYPE'      
      call geti_option(param,config,100,paravi,test)
      if (test.eq.2) then
          do i = 1, test 
             sed_lev(i)=paravi(i)
          enddo  
      else
          call err_option(param,1)  
      endif         
c
      param='-MAG_TYPE'
      call getf_option(param,config,100,paravr,test)
      nlf=test
      if (test.eq.1) then
          mag_lev =  paravr(1)
      else
          call err_option(param,2) 
      endif   
c
      param='-ZMAX_TYPE'
      call getf_option(param,config,100,paravr,test)
      nlf=test
      if (test.eq.1) then
          zmax_lev =  paravr(1)
      else
          call err_option(param,2) 
      endif   
c
      param='-FILT_TYPE'
      call geti_option(param,config,100,paravi,test)
      if (test.ne.1) then
          call err_option(param,3)
      else
         filt_lev=paravi(1)
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
c
      write(*,*) ' filter=',filt_lev
      write(*,*) ' mag=   ',mag_lev
      write(*,*) ' sed=   ',sed_lev(1),sed_lev(2)
      write(*,*) ' Z max= ',Zmax_lev
      write(*,*) ' Input= ',inf(1:lnblnk(inf))
      write(*,*) ' Errors=',errf(1:lnblnk(errf))
      write(*,*) ' Output=',simf(1:lnblnk(simf))
c
cccccccccccc
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  READING THE PHOTOMETRIC LIBRARY and put magnitudes in AB for LF  
c
      write(6,'(A)') '    >>> reading galaxy library  ...',
     > inf(1:(lnblnk(inf)))
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

ccccccccccccccccccccccccccccccccccccccccccccccccc
      open(1,file=simf)
c
cccccccccccccccccccc
      do k = 1, kmax        ! info from input library
c  convert magnitudes from library in AB if needed  
         do j = 1 , nfilt
            if (mag(j,k).ge.999) then 
               mag(j,k)=99.
            else
               if (type(1:1).eq.'V' .or.type(1:1).eq.'v') then 
                  mag(j,k)= mag(j,k) + abcor(j) 
               endif   
            endif
         enddo
c  index range for sed  and redshift
         if (model(k).ge.sed_lev(1).and.model(k).le.sed_lev(2)
     >       .and. z(k).ge.0 .and. z(k).le.zmax_lev ) then
c  compute the color gradients 
           dztot=0
           scale=mag_lev-mag(filt_lev,k)
           do j = 1, nfilt-1
             mag_temp=mag(j,k)+scale
             err_temp= sig_noise(mag_temp,j,Ae,Be,Ce,De)
c             write(*,*) j,mag_temp,err_temp
             do j2=j+1,nfilt
c
             if (mag(j,k).le.90 .and. mag(j2,k).le.90) then 
                dzcol1 =(mag(j,k)-mag(j2,k))
                dzcol2 =(mag(j,k+1)-mag(j2,k+1))
                dzcol =(mag(j,k)-mag(j2,k))
                dzcol = dzcol -(mag(j,k+1)-mag(j2,k+1))
                dzcol = dzcol /dz

                mag_temp=mag(j2,k)+scale
                err_temp2= sig_noise(mag_temp,j2,Ae,Be,Ce,De)
                ecol = sqrt(err_temp**2+err_temp2**2)
c
                dztot=dztot+(dzcol/ecol)**2
                if (model(k).eq.4 .and. z(k).gt.1.9
     >              .and.z(k).le.1.96) then 
c                  write(*,'(2(f5.2,1x),2(I2,1x),7(f10.6,1x))')
c     > z(k),dz,j,j2,scale,mag_temp,dzcol1,dzcol2,dzcol,ecol,dztot
                endif
             endif
             enddo
           enddo  
           dztot=1/sqrt(dztot)
c           if (model(k).eq.4 .and. z(k).ge.1.93
c     >              .and.z(k).le.2.07) then 
c           write(*,*) '---',k,model(k),z(k),dztot
c           endif
           write(1,*) model(k),z(k),dztot
        endif
      enddo  
      close(1)

      STOP
      END


cccccccccccccccccccccccccccccccccc
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
      errmin=0.02
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
      endif  
      sig_noise=dsqrt(sig_noise**2+errmin**2)
c      sig_noise=1.0857/sig_noise 
c
      return
      end         




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











