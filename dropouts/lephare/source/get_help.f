c     last modif 21/09/00
c     AUTHOR : S. ARNOUTS 
c     This subroutine checks if help asked from
c     online option 
c
      subroutine get_help(test)
c
      implicit none 
      integer*4      i,n,iargc,lnblnk,test
      character*4096  var,argv
c
      n= iargc()
      test = 0
c 
      if (n.gt.0) then
        do i = 1, iargc() 
          call getarg(i,argv)
          var=argv
          if (var(1:lnblnk(var)) .eq. '-h' .or. 
     >        var(1:lnblnk(var)) .eq. '-help') then
             test = 1
          endif
        enddo
      endif    
c   
      RETURN
      END
