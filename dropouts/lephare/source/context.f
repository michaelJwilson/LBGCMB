cc++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      FUNCTION bdincl(n,cont,max)
c
      INTEGER FUNCTION  bdincl(n,cont,max)
      implicit none
c 
      integer*4 n,i,max
      real*8    cont,sum    
c
      bdincl=1
      sum = cont
      do i = max,0,-1
         if ((sum-2.**(dble(i))).ge.0. .and. i.gt.n) then
            sum = sum - 2.**dble(i) 
         elseif ((sum-2.**(dble(i))).ge.0. .and. i.eq.n) then
                bdincl=1
                goto 1
         elseif ((sum-2.**(dble(i))).lt.0. .and. i.eq.n) then
                bdincl=0
                goto 1
         endif 
      enddo
c
 1    RETURN
      END
c
