cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     last modif 10/02/2001
c     Author      : S. ARNOUTS 
c     Description : Parabolic interpolation for minmum Chi2
c     Origine     : Bevington book
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c      PROGRAM  int_parab      
      SUBROUTINE int_parab(chxy,chisize,nmax,xb,xmin)
      INTEGER*4 i,nmax,ib,chisize
      REAL*8    chxy(2,chisize),xb
      REAL*8    dx,xmin
c
c search the index value corresponding to z best 
      ib=0
      do i = 1, nmax
         if (DABS(chxy(1,i)-xb).le.1.e-5) then
            ib = i
            goto 10
         endif
      enddo
 10   if (ib.eq.1 .or. ib.eq.nmax) then
         xmin=xb
         goto 20
      endif
      do i = ib-1, ib+1, 1
         if ( chxy(2,i).ge.1e7) then
           xmin=xb
           goto 20
         endif
      enddo   
c compute the best value from the parabolic interpolation of chi2   
      dx=chxy(1,ib)-chxy(1,ib-1)
      if (chxy(2,ib+1)-2*chxy(2,ib)+chxy(2,ib-1).eq.0.d0) then
         xmin =  xb
      else
         xmin=chxy(1,ib+1)-dx*( (chxy(2,ib+1)-chxy(2,ib)) /
     >     (chxy(2,ib+1)-2*chxy(2,ib)+chxy(2,ib-1)) +0.5)
      endif 
c
 20   return
      END
c
