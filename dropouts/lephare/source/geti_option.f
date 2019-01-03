c     last modif 21/09/00
c     description : Return integer value from on-line option 
c           or from config file 
ccccccccccccccccccccccccccccccc
      subroutine geti_option(name,config,nopt,option,test)
cccccccccccccccccccccccccccccc
      implicit none
      integer*4      i,n,iargc,k,lnblnk,test,test0,nopt
      character*4096  var,argv
      character*4096  config,name,name2
      character*4096  paravc(500)
      integer*4       option(500)
c  Option from config file
      test0 = 0
      name2=name(2:lnblnk(name))
      call read_para(name2,config,paravc,test0)
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
            if (nopt.lt.500.and.test.ne.nopt) call err_option(name,3)
         endif   
      enddo
      if (test .eq. 0) test=test0
      do i = 1, 500
         if (i.le.test) then
            read(paravc(i),'(i10)') option(i) 
         else   
            option(i)=0
         endif   
      enddo
c
      RETURN
      END
