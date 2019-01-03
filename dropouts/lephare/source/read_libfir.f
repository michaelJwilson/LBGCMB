c
      subroutine read_libfir(lib,nlib,libfir,size,
     >nfir,liblfir,modlfir,zlfir,lumlfir,reclfir,maglfir,klfir,
     >zso,zmo,dzo,ho,omo,lbdo,mtyp,imag)
c
c      program read_lib
      implicit none 
      integer*4     lnblnk,i,j,k,nrec,size,test,imag,imag1
      integer*4     vali,vali2,recmax,reclibmax,nlib,nbf
      character*4096 zpdir,zpwork,file,valc
      INCLUDE 'dim_filt.decl'      
      INCLUDE 'out_unit.decl' 
      character*4096 lib(5),libfir(5),typ,mtyp
      character*4096 paravc(500),param,emlines
c     FIR library value 
      real*8        zso,zmo,dzo,ho,omo,lbdo
      real*8        zs1,zm1,dz1,h1,om1,lbd1
      integer*4     nfir
      integer*4     modlfir(size),reclfir(size),liblfir(size)
      real*8        zlfir(size),lumlfir(size),maglfir(nbf,size)
      real*8        klfir(nbf,size),val,val2
c      real*8        ebvlfir(size)


c      If you want to give a name to the output screen file
      if(UO.eq.30)
     .      open(UO,file='/tmp/screenLibFir.dat',status='unknown') 


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
c READING LIBRARIES
      nfir=0 
      nrec=0
      do j = 1,nlib  
          file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >            // lib(j)(1:lnblnk(lib(j))) // '.doc'
          param='LIB_TYPE'
          call read_para2(param,file,paravc,test)
          valc=paravc(1)
          if (valc(1:3).ne.'GAL') then 
             write(UO,*) '  FIR library must be GALAXY library'
             nlib=0 
             nfir=0 
             goto 2
          endif
c
          param='NUMBER_ROWS'
          call read_para(param,file,paravc,test)
          read(paravc(1),'(i8)') recmax
          param='RECORD_LENGTH'
          call read_para(param,file,paravc,test)
          read(paravc(1),'(i8)') reclibmax 
c            write(UO,*) recmax,reclibmax
c
          param='FILTERS'
          call read_para2(param,file,paravc,test)
          imag1=test
          if (imag1 .ne. imag) then
             write(UO,*) ' Not right number of filters in FIR LIB'
             write(UO,*) ' -->  Not used '
             nlib=0
             nfir=0  
             goto 2 
          endif  
c
          param='MAG_TYPE'
          call read_para2(param,file,paravc,test)
          typ=paravc(1)
          if ( typ(1:1) .ne. mtyp(1:1) ) then
            write(UO,*) ' Not same  mag system in FIR LIB'
            write(UO,*) ' -->  Not used '
            nlib=0
            nfir=0  
            goto 2  
          endif  
c
          param='LIB_NAME'
          call read_para2(param,file,paravc,test)
          libfir(j) = paravc(1) 
c
          param='Z_STEP'
          call read_para2(param,file,paravc,test)
          read(paravc(1),'(f8.3)') zs1
          read(paravc(2),'(f8.3)') zm1
          read(paravc(3),'(f8.3)') dz1
          if (abs(zs1/zso-1.).gt.0.005 ) then
            write(UO,*) ' Not same Z step  in FIR LIB'
            write(UO,*) ' -->  Not used '
            nlib=0 
            nfir=0 
            goto 2 
          endif  
c
          param='COSMOLOGY'
          call read_para2(param,file,paravc,test)
          read(paravc(1),'(f8.3)') h1
          read(paravc(2),'(f8.3)') om1
          read(paravc(3),'(f8.3)') lbd1
          if (abs(h1/ho-1.).gt.0.02 .OR. ABS(om1/omo-1.).gt.0.02
     >       .OR. ABS(lbd1/lbdo-1.).gt.0.02 ) then
            write(UO,*) ' Not same cosmology  in FIR LIB'
            write(UO,*) ' -->  Not used '
            nlib=0
            nfir=0
            goto 2  
          endif  
c
          param='EM_LINES'
          call read_para(param,file,paravc,test)
          emlines = paravc(1) 
          if (emlines(1:1).eq."Y") then 
              reclibmax =  8*(8+imag1*3)
          else
              reclibmax =  8*(8+imag1*2)
          endif     


          write(UO,*) '     >> FIR Library: ',lib(j)(1:lnblnk(lib(j))),
     >                ' with size:',recmax,reclibmax
          nfir=nfir+recmax
          if (nfir.gt.size) then
             write(UO,*) ' dim lib :',nfir,'must be < :',size
             STOP
          endif   
c reading the library 
          file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >           // lib(j)(1:lnblnk(lib(j))) // '.bin'
          open(2,file=file,status='unknown',access='direct',
     >               recl=reclibmax)
c          write(UO,*) 'reading file ...'
          do i = 1, recmax 
             nrec=nrec+1
             liblfir(nrec)=j     
c             write(UO,*) i,nrec 
             if (emlines(1:1) .eq. "Y") then   
                read(2,rec=i) modlfir(nrec),vali2,val,lumlfir(nrec),
     >                      zlfir(nrec),val,reclfir(nrec),
     >    vali,(maglfir(k,nrec),k=1,vali),(klfir(k,nrec),k=1,vali),
     >                                (val2,k=1,vali)
             else
                read(2,rec=i) modlfir(nrec),vali2,val,lumlfir(nrec),
     >                      zlfir(nrec),val,reclfir(nrec),
     >    vali,(maglfir(k,nrec),k=1,vali),(klfir(k,nrec),k=1,vali)
             endif
          enddo
c 1        close(2) 
          close(2) 
c
      enddo
c

c
 2    continue 
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      RETURN
 56   write (UO,*) 'File ',file(1:lnblnk(file)),' not found -> STOP '
      stop
      END
