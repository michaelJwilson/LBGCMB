c      last modif 28/10/2000
c      Author : S. ARNOUTS 
c      description : Extract values from character string
c      if not found returns test=0
ccccccccccccccccccccccccccccccccccccccc 
      subroutine get_value(string,paravc,test)
      implicit none 
c     
      integer*4     i,imin,imax,ilen,inum,test
      character*4096  string
      character*4096  var
      character*4096  paravc(500)
c
c     search the first non blank character and the last character
      imin=0
      imax=0
      do i = 1, 500
         paravc(i) = '       ' 
      enddo   
      do i = 1, 4096
         var=string(i:i)
         if (var(1:1).ne.' '.and.var(1:1).ne.char(9)
     >       .and.imin.eq.0) then 
            imin=i
         elseif ((var(1:1).eq.' '.or.var(1:1).eq.'#'.or.
     >        var(1:1).eq.char(9)).and.imin.gt.0) then
            imax=i-1
            goto 1
         endif
      enddo
 1    if (imax.lt.imin) then
         test=0
         return
      endif   
      var=string(imin:imax)
      string=var
c
c extract the different elements of the sub-string
      test=0
      ilen=imax-imin+1
      inum=1
      do i =1 , ilen
         var=string(i:i)
         if (var(1:1).eq.','.and.i.ne.ilen) then
            test=test+1
            paravc(test)=string(inum:i-1)
            inum=i+1
         endif
      enddo   
      test=test+1
      paravc(test)=string(inum:ilen)
c
      return
      end
