c       
        PROGRAM test
 
       implicit none 
       integer*4  k,i
       integer*4  nr,nmod,iw
       real*8     dummy,vec1,vec2,vec     
       character*256  file


       file='/Users/arnouts/zpwork_dev/lib_bin/LIB_BCSTOCH.bin'
       file=file(1:lnblnk(file))
c
        open(1,file=file,status='unknown',access='direct',
     >                                  recl=55368)
c
        read(1,rec=1) nr,nmod,dummy,iw,(vec,k=1,iw)
            write(*,*) nr,nmod,dummy,iw,vec
        read(1,rec=2) nr,nmod,dummy,iw,vec1,vec2,(vec,k=3,iw)
            write(*,*) nr,nmod,dummy,iw,vec1,vec2,vec
       read(1,rec=3) nr,nmod,dummy,iw,(vec,k=1,iw)
            write(*,*) nr,nmod,dummy,iw
       read(1,rec=4) nr,nmod,dummy,iw,(vec,k=1,iw)
            write(*,*) nr,nmod,dummy,iw
       read(1,rec=5) nr,nmod,dummy,iw,(vec,k=1,iw)
            write(*,*) nr,nmod,dummy,iw
       read(1,rec=6) nr,nmod,dummy,iw,(vec,k=1,iw)
            write(*,*) nr,nmod,dummy,iw
       read(1,rec=7) nr,nmod,dummy,iw,(vec,k=1,iw)
       do i = 1,100002
         read(1,rec=i) nr,nmod,dummy,iw,(vec,k=1,iw)
         write(*,*) nr,nmod,dummy,iw,vec
       enddo        
       close(1)
c
       end
