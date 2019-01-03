c   last modif : 23/07/002
      program read_specz0
c   
c     this program extracts the spectrum for Galaxy 
c     from the library and put it at the obs. z and z=0 and used extinction ...
c     FOR HDF N + S  in order to run PCA analysis  
      implicit none
      integer*4  imas,i,j,k,iw,modb,iext,km
      integer*4  wmax,test,model
      integer*4  recmax,nrec(10000),reclmax,recmax0,recmin
      integer*4  lnblnk,nmod,nop,nr,pass,sampl,iopa(81)
c      
c      parameter  (wmax=8000)
      INCLUDE 'dim_wave.decl'      
c
      real*8     wave(wmax),flz0(wmax,1000)
      real*8     exti(2,wmax),opal(81,wmax),opat(81,wmax),extis(2,wmax)
      real*8     dummy,vec(wmax),ebv,wav(wmax),fl(wmax)
      real*8     c,ageb,zb,val,chib
      parameter(c=2.99792458e+18) 
      character*512  file,lib,zpdir,extlaw,param,zpwork
      character*512  paravc(100),fileop(100)
      character*1024 str
      recmin=0
      recmax=0
      sampl=1
c  fixing the libraries
c
      extlaw = 'calzetti.dat'
      lib  = 'LIB_GIS'      
c
c  environmental variable 
      call getenv('LEPHAREDIR',zpdir)
      test=lnblnk(zpdir)
      if (test .eq. 0) then
        write(*,*) 'WARNING :  variable LEPHAREDIR not defined'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(6,*) 'WARNING :  variable LEPHAREWORK not defined'
        stop
      endif
c
      do k = 1, wmax 
         wav(k)=0
         fl(k)=0
      enddo
c  read extinction 
c      
        file = zpdir(1:lnblnk(zpdir)) //'/ext/' 
     >   // extlaw(1:lnblnk(extlaw))
        open(1,file=file,status='unknown')
        i=0
        do while (.true.)
          i=i+1
          read(1,*,end=11) exti(1,i),exti(2,i)
        enddo
 11     iext=i-1
        close(1)      
c  read opacity
c
        file = zpdir(1:lnblnk(zpdir)) // '/opa/OPACITY.dat'
        open(1,file=file,status='unknown')
        do i = 1,81
           read(1,'(A)') str 
           call val_string(str,paravc,test)
           fileop(i)=paravc(2)
        enddo
        close(1)
        do i = 1, 81
           file = zpdir(1:lnblnk(zpdir)) // '/opa/' 
     > // fileop(i)
           open(1,file=file,status='unknown')
           k = 0
           do while (.true.)
              k = k + 1
              read(1,*,end=12) opal(i,k), opat(i,k)
           enddo  
 12        iopa(i) = k - 1
           close(1)
        enddo        
c  loop on input file to extract spectrum 
c
      open(35,file='hdfns_spec.dat')   ! output with spectra at z=0
      open(33,file='hdfns_mod.list') 
      read(33,*) 
      read(33,*) 
      km = 0 
      do while (.true.)
         km = km + 1
         read(33,*,end=13) modb,ageb,ebv,chib,zb    
c         write(6,*) km,'--  objects --> ',modb,ageb,ebv,zb
c
c  read number of CSP model and number of records per model nrec(i)
c      write(*,*) lib(1:lnblnk(lib))
      file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     > // lib(1:lnblnk(lib)) // '.doc'

      param='NUMBER_SED'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') imas
      param='RECORD_LENGTH'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') reclmax
      recmax0 = 0
      do i = 1, imas
         write(str,'(i6.6)') i 
         param = 'MOD_' // str(1:lnblnk(str))
         call read_para2(param,file,paravc,test)
         read(paravc(1),'(i8)') model
         if (test.ne.0) then
            read(paravc(3),'(i8)') nrec(i)
           if (modb.eq.model) then
             recmin=recmax0+1
             recmax=recmin+ nrec(i)
             goto 30
           else 
             recmax0 = recmax0 + (nrec(i) + 1)
           endif
         else
            call err_option(param,1)
         endif   
      enddo   
c   ... READING MODELS  !!
 30   file = zpwork(1:lnblnk(zpwork)) // '/lib_bin/'
     >       // lib(1:lnblnk(lib)) // '.bin'
      open(28,file=file,status='unknown',access='direct',recl=reclmax)
      pass=0
      val = zb/0.1 + 1
      nop = int(val)
c
      do j = recmin, recmax
         read(28,rec=j) nr,nmod,dummy,iw,(vec(k),k=1,iw)
         if (dummy.eq.-1.0) then
c              write(*,*) km,nmod,' lambda ...'
              do k = 1,iw
                 wave(k) = vec(k)
              enddo
         else
            if (DABS((dummy-ageb)/1e9).le.0.001) then
              write(6,*) km,'-- spectrum out <--',modb,dummy, ageb
c               call scale_opa(wave,vec,iw,opal,opat,iopa,zb,wn,fn,inmax)
              if (ebv.gt.0) then
                  call lambda(wave,iw,exti,iext,extis)
              endif   
c               do k = 1, inmax ,sampl
              do k = 1, iw ,sampl
c  at z = 0 
                  flz0(k,km) = vec(k)
                  if (ebv.gt.0) then
                     flz0(k,km)=flz0(k,km)*10**(-0.4*ebv*extis(2,k))
                  endif   
c at z = zb
c                  wav(k)=wave(k)*(1+zb)
c                  fn(k)=vec(k)/(1+zb)
c                  if (ebv.gt.0) then
c                     fn(k)=fn(k)*10**(-0.4*ebv*extis(2,k))
c                     write(*,*) ebv, extis(2,k),fn(k)
c                  endif   
              enddo 
              do k = 1, iw
                write(35,'(I3,1x,E12.6,1x,E12.6)') km,wave(k),flz0(k,km)
              enddo
c              write(6,*) nmod,dummy,iw,' flux...'
              pass=1
              goto 32 
            elseif (dummy.ne.ageb .and. pass.eq.1) then 
               goto 32
            endif   ! age
         endif      ! dummy
      enddo         ! record
      write(6,*) ' WARNING : No spectrum found ...'
 32   close(28)
c
      enddo         ! list 
 13   close(33)
      close(35)
      km = km - 1
c format for PCA 
      open(1,file='hdf_pca.dat')
      write(1,*) ((flz0(k,i),k=1,iw),i=1,km)
      close(1)
      write(6,*) 'number of spectra :',km
      write(6,*) 'number of wavelength :',iw
      write(6,*) '    '
      write(6,*) 'spectra in output file hdfns_spec.dat '
      write(6,*) 'spectra PCA in  hdf_pca.dat  '
      write(6,*) '    '
c
      STOP
      end
c
