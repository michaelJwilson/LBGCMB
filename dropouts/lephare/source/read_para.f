c      last modif 06/01/2000
c      Author : S. ARNOUTS 
c      description : Read input and output parameter "name" for zphot
c      and extract its value.  If not found returns test=0
ccccccccccccccccccccccccccccccccccccccc 
      subroutine read_para(name,conf,paravc,test)
      implicit none 
c     
      integer*4      i,j,imax,imin,lnblnk,test,index
      character*4096  name
      character*4096  var,string
      character*4096  var2
      character*4096  conf
      character*4096  paravc(500)
c
      test=0
      open(1,file=conf,status='old',err=56)
      i=0
      imax=1
      do while(.true.)
         i = i + 1
         read(1,'(a)',end=10) var
         do j=1,4096
            var2=var(j:j)
            if (var2(1:1).eq.' ' .or. var2(1:1).eq.char(9)) then
               imax=j-1
               goto 5
            endif   
         enddo
 5       if (var(1:imax).eq.name(1:lnblnk(name))) then
           test=1
c extract name to var 
            imin=index(var,'name(1:lnblnk(name))')
            imax=lnblnk(name)+imin+1
            string=var(imax:lnblnk(var))
c get the value for parameters
            call get_value(string,paravc,test)           
           goto 10 
         endif
      enddo
 10   close(1)


      return

 56   write (6,*) 'File ',conf(1:lnblnk(conf)),' not found -> STOP '
      stop
c
      end
