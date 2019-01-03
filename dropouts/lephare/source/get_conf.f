c     last modif 21/09/00
c     AUTHOR : S. ARNOUTS 
c     This subroutine extracts the config file to be used from 
c     online option 
c
      subroutine get_conf(name,out,test)
c
      implicit none 
      integer*4      i,n,iargc,j,lnblnk,test
      character*4096  name,out,var,argv
c
      n= iargc()
      test = 0
      if (n.gt.0) then
        do i = 1, iargc() 
          call getarg(i,argv)
          var=argv
          if (var(1:lnblnk(var)) .eq. name(1:lnblnk(name))) then
             j=i+1
             call getarg(j,argv)
             out=argv
             test = 1
          endif
        enddo
      endif    
c   
      RETURN
      END
