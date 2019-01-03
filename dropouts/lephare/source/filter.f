      program FILTER
c     Create a unique filter file (FILTER_FILE) with 
c     all the single files in the list FILTER_LIST
c 
      implicit none
      integer*4      i,j,nbf,test,lnblnk,nfile,imax
      integer*4      index,idot,maxsize
      INCLUDE 'out_unit.decl' 
      INCLUDE 'dim_filt.decl'
      parameter (maxsize=110000)
      real*8         lbd,trans,transf,transtyp(500),paravr(500),lbmean
      real*8         lb(maxsize),tr(maxsize),area
      integer*4       calibtyp(500),paravi(500)
      character*4096  file(nbf),filters,zpdir,config,filtfile
      character*4096  param,filtdoc,filt,zpwork
      character*4096  paravc(500),name
      character*4096 str
      real*8     zp(nbf),flmoy(nbf),flwidth(nbf),abcor(nbf)
      real*8     fcorr(nbf)

c      If you want to give a name to the output screen file
      if(UO.eq.30)
     .      open(UO,file='/tmp/screenSedtolib.dat',status='unknown') 

      
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
      param='filter'
      call get_help(test)
      if (test .eq. 1) call help(param) 
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      param='-c'
      call get_conf(param,config,test)
      if (test.ne.1) call err_option(param,1)
      call get_path(config)
c      
c  read option
c
      param='-FILTER_LIST'
      call getc_option(param,config,500,paravc,test)
      if (test.ge.1) then 
         do i = 1,test
           file(i)=paravc(i)(1:lnblnk(paravc(i)))
         enddo  
         nfile=test
      else
         call err_option(param,1)
      endif
      if (nfile.gt.nbf) then 
        write(UO,*) 'Number of filters ',nfile,' larger than the number'
        write(UO,*) ' allocated in  dim_filt.decl :',nbf
        write(UO,*) ' increase values in dim_filt.dec --> STOP ' 
        stop
      endif  
c
      param='-TRANS_TYPE'
      call getf_option(param,config,500,paravr,test)
c      write(*,*) test
      if (test.eq.0) then
         do i = 1, nfile 
            transtyp(i) = 0
         enddo
         call err_option(param,4)
      elseif (test.eq.1) then
         if (paravr(1).ne.0. .and. paravr(1).ne.1.) paravr(1)=0. 
           do i = 1, nfile
              transtyp(i)=paravr(1)
           enddo   
           
      elseif (test.eq.nfile) then
         do i = 1, nfile
         if (paravr(i).ne.0. .and. paravr(i).ne.1.) paravr(i)=0. 
            transtyp(i)=paravr(i)
         enddo   
      elseif (test.ne.nfile) then
         do i = 1, nfile 
            transtyp(i) = 0
         enddo
         call err_option(param,4)         
      endif  
c
      param='-FILTER_CALIB'
      call geti_option(param,config,500,paravi,test)
      if (test.eq.0) then
         do i = 1, nfile 
            calibtyp(i) = 0
         enddo
         call err_option(param,4)
      elseif (test.eq.1) then
c         if (paravi(1).ne.0 .and. paravi(1).ne.1) paravi(1)=0 
           do i = 1, nfile
              calibtyp(i)=paravi(1)
           enddo             
      elseif (test.eq.nfile) then
         do i = 1, nfile
c         if (paravi(i).ne.0 .and. paravi(i).ne.1) paravi(i)=0 
            calibtyp(i)=paravi(i)
         enddo   
      elseif (test.gt.1 .and. test.ne.nfile) then
         do i = 1, nfile 
            calibtyp(i) = 0
         enddo
         call err_option(param,4)         
      endif  

c
      param='-FILTER_FILE'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  filters=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
      idot=index(filters,'.')
c
      filtfile=zpwork(1:lnblnk(zpwork)) //'/filt/' // 
     > filters(1:lnblnk(filters))
c      
      filtdoc=zpwork(1:lnblnk(zpwork)) // '/filt/' //    
     > filters(1:idot-1) // '.doc'
c
c   INFO ON SCREEN
        write(UO,'(A)') "#####################################"
        write(UO,'(A)') "# It s building the filter file     #"
        write(UO,'(A)') "#   with the following options :    #"
        write(UO,'(2A)')"# Config file : ",config(1:lnblnk(config))
        write(UO,'(A,500(1x,f2.0))') "# TRANS_TYPE  : ",
     >                                (transtyp(i),i=1,nfile)
        write(UO,'(A,500(1x,I2))') "# FILTER_CALIB: ",
     >                                (calibtyp(i),i=1,nfile)

        write(UO,'(2A)')"# FILTER_FILE : ",filtfile(1:lnblnk(filtfile))
      write(UO,'(2A)')"# FILTER_FILE.doc: ",filtdoc(1:lnblnk(filtdoc))
        write(UO,'(500A)')"# FILTER_LIST : ",(file(i)
     >               (1:lnblnk(file(i))),' ',i=1,nfile)
        write(UO,'(A)') "#####################################"
c
      open(1,file=filtfile,status='unknown')
      open(2,file=filtdoc,status='unknown')
      write(1,100) nfile
      do i = 1,nfile
         filt=zpdir(1:lnblnk(zpdir)) //'/filt/'//
     >   file(i)(1:lnblnk(file(i)))
c         write(*,*) i,nfile,filt
         open(3,file=filt,status='old',FORM='formatted',err=56)
         read(3,'(A)') str
         call val_string(str,paravc,test)
c  extract the name of the filter ...
         name=paravc(2)
c  get the number of row for the filter 
         j=0
         lbmean = 0
         do while (.true.) 
          j=j+1
          read(3,*,end=10)  lb(j), tr(j)
         enddo
 10      imax=j-1  
         close(3)
         area = 0
         lbmean=0  
         do j = 1,imax-1
           trans =(tr(j)+tr(j+1))/2.
           area  =area + trans*(lb(j+1)-lb(j))
           lbmean=lbmean+trans*(lb(j+1)-lb(j))*(lb(j)+lb(j+1))/2.
         enddo
         lbmean = lbmean / area
c put nb row , name and transmission in final filter file 
         write(1,110) imax,name,calibtyp(i),i
         write(2,110) imax,name,calibtyp(i),i
         write(*,110) imax,name(1:lnblnk(name)),calibtyp(i),i
         open(3,file=filt,status='old',FORM='formatted')
         read(3,'(A)') str
         do j =  1, imax
           read(3,*) lbd,trans
           transf = trans*(lbd/lbmean)**transtyp(i)
           write(1,*) lbd,transf,i
         enddo
         close(3)
      enddo
      close(1)
      close(2)
      if(UO.eq.30) close(30)
c
      call zeropoint(filters,zp,abcor,flmoy,flwidth,fcorr,nfile)
c
 100  format("# ",1x,I5)
 110  format("# ",1x,I5,1x,A15,1x,I5,3x,I2)

      stop
 56   write (6,*) 'File ',filt(1:lnblnk(filt)),' not found -> STOP '
      
      end
