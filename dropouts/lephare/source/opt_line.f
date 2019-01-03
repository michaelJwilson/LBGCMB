c
cccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine opt_line(name,nopt,option,test)
      implicit none
      integer*4     i,n,iargc,j,k,lnblnk,test,test0,nopt
      character*4096 var,argv
      character*4096 name
      character*4096 option(500),paravc(500)
c  Option from command line 
      test0=0
      test=0
      n= iargc()
      do i = 1, n
         call getarg(i,argv)
         var=argv
         if (var(1:lnblnk(var)) .eq. name(1:lnblnk(name))) then
            k = i + 1
            call getarg(k,argv)
            var=argv
            call get_value(var,paravc,test)
            if (nopt.lt.500 .and. test.ne.nopt) then
               call err_option(name,3)
            else
               do j = 1,nopt
                  option(j) = paravc(j)
               enddo
            endif   
         endif   
      enddo
      if (test .eq. 0) then 
        test=test0
      else
        do i = test+1,500
           option(i)=' '
        enddo
      endif
      RETURN
      END
