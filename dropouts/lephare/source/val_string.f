c      last modif 28/10/2000
c      Author : S. ARNOUTS 
c      description : Extract values from character string
c      if not found returns test=0
ccccccccccccccccccccccccccccccccccccccc 
      subroutine val_string(string,paravc,test)
      implicit none 
c     
      integer*4      i,imin,imax,inum,test,lnblnk,bl(4096)
      character*4096  string
      character*4096  var
      character*4096  paravc(500)
c
c     search the first non blank character and the last character
      imin=0
      imax=0
      do i = 1, 500
         paravc(i) = '            ' 
      enddo   
c
c     extract the different elements of the sub-string
      do i = 1, 4096
         bl(i)=0
      enddo   
      do i =1 , lnblnk(string)
         var=string(i:i)
         if (var(1:1).eq.','.or.var(1:1).eq.' '.or.
     >        var(1:1).eq.char(9) ) then
             bl(i)=i
         endif
      enddo
      test=0
      inum=0
      if (bl(1).ne.0) then
         inum=1
      endif   
      bl(lnblnk(string)+1)=lnblnk(string)+1
      do i = 2, lnblnk(string)+1
         if (bl(i).ne.0 .and. bl(i-1).eq.0) then
            if (test.lt.499) then
               test=test+1
               paravc(test)=string(inum+1:i-1)
               inum=i
            elseif (test.eq.499) then
               test=test+1
               paravc(test)=string(inum+1:lnblnk(string)+1)
               return
            else
c              write(6,*) "WARNING: Keyword with more than 500 values"
            endif 
         elseif (bl(i).ne.0) then
            inum=i
         endif
      enddo   
c
      return
      end















