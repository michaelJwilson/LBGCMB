c
       subroutine  check_libphys(libphys,nlibphys,libname,libppara,
     >      zrecp,zrecpi,zrecps,nzrecp,
     >      h0,om0,l0,mtyp,imag) 

      implicit none 
      integer*4      maxsize
      parameter     (maxsize=110000)
      INCLUDE 'out_unit.decl' 
      integer*4      i,j,imax,lnblnk
      integer*4      k,kmax,kphys,recp,modp
      character*4096  param,valc,var,var2,file,mtyp
      character*4096  paravc(500),zpdir,zpwork,libphys,libname
      integer*4       imag,test,nlibphys,isel
      real*8          h0,om0,l0,hub,omeg,lbda
      real*8          libppara(50,maxsize),zrecp(500)
      integer*4       zrecpi(500),zrecps(500),nzrecp
c    
        if(UO.eq.30)
     >      open(UO,file='/tmp/screenReadLib.dat',status='unknown') 
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
cccccc    reading doc file
        file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >            // libphys(1:lnblnk(libphys)) // '.doc'
        write(UO,*) '   checking the Physical library :',
     >    file(1:lnblnk(file))
cccccc
        param='LIB_NAME'
        call read_para2(param,file,paravc,test)
        libname=paravc(1)
        param='MAG_TYPE'
        call read_para2(param,file,paravc,test)
        valc=paravc(1)
        if (valc(1:lnblnk(valc)).ne.mtyp(1:lnblnk(mtyp))) then 
          write(UO,*) ' WARNING: MAG_TYPE differs between library'
          write(UO,*) '  -->   PHYSICAL Library switches OFF'
          nlibphys=0
        endif
        param='FILTERS'
        call read_para2(param,file,paravc,test)
        if (test.ne.imag) then 
          write(UO,*) ' WARNING: nb of filters differs between library'
          write(UO,*) '  -->   PHYSICAL Library switches OFF'
          nlibphys=0
        endif
        param='COSMOLOGY'
        call read_para2(param,file,paravc,test)
        read(paravc(1),'(f8.3)') hub
        read(paravc(2),'(f8.3)') omeg
        read(paravc(3),'(f8.3)') lbda
        if (DABS(h0/hub-1.)   .gt. 0.02  .or. 
     >      DABS(om0/omeg-1.) .gt. 0.02  .or.
     >      DABS(l0/lbda-1.)  .gt. 0.02       ) then   
          write(UO,*) ' WARNING: Cosmology  differs between library'
          write(UO,*) '  -->   PHYSICAL Library switches OFF'
          nlibphys=0
        endif
cccccccccccccccccccccccccccccccccccccccccccc
c   read the redshift record 
      open(1,file=file(1:lnblnk(file)),status='old',err=56)
      test=0
      isel=0
      imax=1
      i=0
      do while(.true.)
         i = i + 1
         read(1,'(a)',end=10) var
         do j=1,4096
            var2=var(j:j)
            if (var2(1:1).eq.' ' .or. var2(1:1).eq.char(9)) then
               imax=j-1
               goto 1
            endif   
         enddo  
 1       if (var(1:imax).eq."Z_RECORD") then
           test=1
           isel= i
           goto 10
         endif
      enddo
 10   close(1)

      if (test.eq.0 .or. isel.eq.0) then 
          write(UO,*) ' WARNING: Z_RECORD not found in PHYS library'
          write(UO,*) '  -->   PHYSICAL Library switches OFF'
          nlibphys=0
      else
         open(1,file=file(1:lnblnk(file)),status='old',err=56)
         do i = 1,isel-1
            read(1,*)
         enddo
         read(1,*) var,nzrecp,(zrecp(i),zrecpi(i),zrecps(i),i=1,nzrecp)
         close(1)
      endif 

cccccccccccccccccccccccccccccccccccccccccccc
cccccc    reading phys file
       write(UO,*) ' extracting physical parameters ...'
       file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >        // libname(1:lnblnk(libname)) // '.doc'
c
       param='NUMBER_ROWS'
       call read_para(param,file,paravc,test)
       read(paravc(1),'(i8)') kmax
c
       if (kmax .gt. 0 ) then 
          file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >        // libname(1:lnblnk(libname)) // '.phys'
          open(1,file=file,status='unknown',err=56)
          k = 0
          do while (.true.)
             k = k + 1
             read(1,*,end=12) recp,modp,(libppara(i,k),i=1,31)
          enddo  
 12       kphys = k -1
          close(1)
          write(UO,*) 'number of physical parameters=',kphys,kmax
          if (kphys.ne.kmax) then 
         write(UO,*) ' Size problem for the  lib_bin library ',
     >                 libname(1:lnblnk(libname)),kphys,kmax          
            write(*,*) ' check the lib_bin *.phys files ...'
            stop
          endif
c
          write(*,*) "1 , 0 ",libppara(1,1),libppara(6,1),libppara(31,1)
          write(*,*) "2 , 1 ",libppara(1,2),libppara(6,2),libppara(31,2)
          write(*,*) "3 , 2 ",libppara(1,3),libppara(6,3),libppara(31,3)
          write(*,*) "10000, 9999 ",
     >     libppara(1,10000),libppara(6,10000),libppara(31,10000)
          write(*,*) "10001, 10000 ",
     >     libppara(1,10001),libppara(6,10001),libppara(31,10001)

          write(*,*) "100000, 99999 ",
     >     libppara(1,100000),libppara(6,100000),libppara(31,100000)
          write(*,*) "100001, 100000 ",
     >     libppara(1,100001),libppara(6,100001),libppara(31,100001)
       endif 
c
c
        return
 56   write (UO,*) 'File ',file(1:lnblnk(file)),' not found -> STOP '
      stop

        end
