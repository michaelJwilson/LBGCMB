c     last modif 21/09/00
c     description : Return Character value from on-line option 
c           or from config file 
ccccccccccccccccccccccccccccccc
      subroutine getc_option(name,config,nopt,option,test)
cccccccccccccccccccccccccccccc
      implicit none
      integer*4      i,n,iargc,j,k,lnblnk,test,test0,nopt
      character*4096  var,argv
      character*4096  config,name,name2
      character*4096  option(500),paravc(500)
c  Option from config file
      test0 = 0
      name2=name(2:lnblnk(name))
      call read_para(name2,config,option,test0)
c  Option from command line 
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
