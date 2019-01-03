      program read_libphys
      implicit none
      character*4096 file,a
      integer*4      reclibmax,nrec,i,j,k,iform,irec
      integer*4      nmod,iel,nr,imag,ifilt
      real*8         valebv,zoss,eta,mag(50),kcor(50)
      real*8         zrec(3,1000),z(10000000)
      real*8         zr(1000)
      integer*4      reci(1000),recs(1000)
      reclibmax=168
      nrec=21225
      file='/Users/sarnouts/zpwork_dev/lib_mag/BCSTOCH2_HDF.bin'
      nrec=1560054
      file='/Users/sarnouts/zpwork_dev/lib_mag/BCSTOCH_HDF.bin'
      open(2,file=file,status='unknown',access='direct',
     >               recl=168)

      write(*,*) file(1:lnblnk(file))
      k=0
      do j = 1,nrec
        read(2,rec=j)  nmod,iel,valebv,z(j),
     >                eta,nr,imag,(mag(ifilt),ifilt=1,imag),
     >                (kcor(ifilt),ifilt=1,imag)

        if (j.eq.1) then
           k=k+1 
           zrec(1,k)=z(j)
           zrec(2,k)=j
        elseif (z(j).gt.z(j-1)) then
           zrec(3,k)=j-1
           write(*,'("Z_RECORD  ",f6.3,i12,1x,i12,2x,I12)')
     >   zrec(1,k),(idint(zrec(i,k)),i=2,3),idint(zrec(3,k)-zrec(2,k)+1)
           k=k+1
           zrec(1,k)=z(j)
           zrec(2,k)=j
        endif 
        if (j.eq.nrec) zrec(3,k)=j 
      enddo
      write(*,'("Z_RECORD  ",f6.3,i12,1x,i12,2x,I12)')
     >  zrec(1,k),(idint(zrec(i,k)),i=2,3),idint(zrec(3,k)-zrec(2,k)+1)
      close(2)
      iform=k 

       open (1,file='temp',status='unknown')  
        write(1,'(A,1x,I6,1x,500(f6.3,2x,2(I12,2x)))') 'Z_RECORD ',iform
     >  ,(zrec(1,k),(idint(zrec(j,k)),j=2,3),k=1,iform) 
       close(1)
c
       write(*,*) 'reading temporary file temp  ...'
       open (1,file='temp',status='unknown')  
       read(1,*) a,irec,(zr(k),reci(k),recs(k),k=1,irec) 
       close(1) 
       do k = 1,irec 
         write(*,*)  zr(k),reci(k),recs(k),(recs(k)-reci(k)+1)
       enddo


      stop
      end
