c     last modif 21/09/00
c     description : Return floating value from on-line option 
c           or from config file 
ccccccccccccccccccccccccccccccc
      subroutine getf_option(name,config,nopt,option,test)
cccccccccccccccccccccccccccccc
      implicit none
      integer*4       i,n,iargc,j,k,lnblnk,imin,idot,test,test0,nopt
      character*4096  var,argv
      character*4096  config,name,name2
      character*4096  paravc(500)
      real*8          option(500)
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
            k= i + 1
            call getarg(k,argv)
            var=argv
            call get_value(var,paravc,test)
         if (nopt.lt.500.and.test.ne.nopt) call err_option(name,3)
         endif   
      enddo
      if (test .eq. 0)  test=test0
      do i = 1,500
         if (i.le.test) then
c  check if a "." is detected to avoid problem in float conversion
            imin=0
            idot=0
            do j = 1, lnblnk(paravc(i))+1
               var=paravc(i)(j:j)
               if (var(1:1).eq.'.') then
                  idot=1
               endif   
               if (var(1:1).eq.'e' .or. var(1:1).eq.'E'
     >        .or. var(1:1).eq.'d' .or. var(1:1).eq.'D') imin=j 
            enddo
            if (idot.eq.0 .and.imin.eq.0) then
               paravc(i) = paravc(i)(1:lnblnk(paravc(i))) // '.'
            elseif (idot.eq.0 .and.imin.ne.0) then 
               paravc(i) =  paravc(i)(1:imin-1) 
     >             // '.' // paravc(i)(imin:lnblnk(paravc(i)))
            endif   
            read(paravc(i),'(E23.12E3)') option(i)
         else       
           option(i)=0.
         endif  
      enddo
      RETURN
      END



