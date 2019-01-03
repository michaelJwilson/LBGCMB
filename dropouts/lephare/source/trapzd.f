c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Calcul l'aire sous la fonction y entre x=a et x=b
c
c      PROGRAM TRAPZD
      SUBROUTINE TRAPZD(x,y,bmax,a,b,s)
      implicit none 
c      
      INTEGER*4 bmax,i,k,imin,imax,pass
      REAL*8    a,b,x(bmax),y(bmax)
      REAL*8    xmin,ymin,xmax,ymax
      REAL*8    s,sum

c      a = 8.
c      b = 1000.
c
c   Initialisation
      k    = 0
      imax = 0
      s    = 0.
      sum  = 0.  
      pass = 0
      imin = 1
      imax = bmax
      xmin = x(1)
      xmax = x(bmax)
      ymin = y(1)
      ymax = y(bmax)
c      write(*,*) xmin,xmax,ymin,ymax
c  Check if a and b defined in x(i) step 
      do i = 1, bmax
         if (DABS(x(i)-a).lt.1.e-5) then
            xmin= a
            ymin= y(i)
            imin= i
            pass=pass+ 1
            goto 10
         endif
      enddo   
 10   do i = imin, bmax
         if (DABS(x(i)-b).lt.1.e-5) then
            xmax= b
            ymax= y(i)
            imax= i
            pass=pass+2
            goto 20 
         endif
      enddo   
 20   if (pass.eq.3) goto 40
c  Here Search for the limits if a or b not defined in x(i), y(i)
      if (pass.eq.0 .or. pass.eq.2) then   ! if a not defined in x(i)
        do i = 1, bmax-1
          if (x(i).lt.a .and. x(i+1).gt.a) then
            xmin= a
            ymin= y(i) + (y(i+1)-y(i))/dabs(x(i+1)-x(i))*dabs(xmin-x(i))
            imin= i
            goto 30
          endif
        enddo            
      endif
 30   if (pass.eq.0 .or. pass.eq.1) then ! if b not defined in x(i)
        do i = 1, bmax-1
          if (x(i).lt.b .and. x(i+1).gt.b) then
            xmax= b
            ymax= y(i) + (y(i+1)-y(i))/dabs(x(i+1)-x(i))*dabs(xmax-x(i))
            imax= i
            goto 40
          endif
        enddo      
      endif   
c  Integration of the inner part excluding limits                 
 40   continue
      do  i=imin+1,imax-2
            sum = sum + (x(i+1)-x(i))*(y(i)+y(i+1))/2.
      enddo    
c  Add the area with borders  
        s=sum+ dabs(x(imin+1)-xmin)*(ymin+y(imin+1)) /2.
        s=s  + dabs(xmax-x(imax-1))*(ymax+y(imax-1)) /2.
c
c        write(*,'(3(I4,1x),2(E12.5,1x))') pass,imin,imax,a,b
c        write(*,'(I4,5(E12.5,1x))') pass,xmin,ymin,xmax,ymax,s 
c
      return
      END

c
