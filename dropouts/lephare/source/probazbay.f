c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Return zmin,zmax included inside 50%,68%,90%,99% of the area 
c     zpdzi-s(1) : +/- 50%
c     zpdzi-s(2) : +/- 68%
c     zpdzi-s(3) : +/- 90%
c     zpdzi-s(4) : +/- 99%
c
      SUBROUTINE PROBAZBAY(x,y,bmax,area,zmed,zpdzi,zpdzs)
      implicit none 
c      
      INTEGER*4 bmax,i
      REAL*8    x(bmax),y(bmax),yn(bmax),zmed,slope
      REAL*8    zli,zls,area,proba,probai(bmax),zpdzi(10),zpdzs(10)
c
c   Initialisation
      zli   = x(1)
      zls   = x(bmax)
      proba = 0.
      do i = 1, bmax 
         yn(i) = y(i)/area
         if (yn(i) .le. 1.d-60) yn(i)=0.
      enddo      
c
      zmed= -99
      do i = 1, 10
        zpdzi(i) = -99.
        zpdzs(i) = -99.
      enddo  
c
      probai(1)=0.d0
      do i = 1, bmax-1
         zls=x(i+1)
         call TRAPZD(x,yn,bmax,zli,zls,proba)
         probai(i+1) = proba
c         write(*,*) zli,zls,probai(i)
      enddo
      do i = 1, bmax-1
c   Zmedian 
         if (probai(i).lt.0.5 .AND. probai(i+1).ge.0.5) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zmed=x(i) + (0.5d0 - probai(i))/slope
c            if (DABS(probai(i+1)-0.5d0).lt.0.00001d0) then
c               zmed = x(i+1) 
c            else
c               zmed = (x(i+1)+x(i))/2.d0 
c            endif
         endif         
c   50% inf  
         if (probai(i).lt.0.25 .AND. probai(i+1).ge.0.25) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zpdzi(1)=x(i) + (0.25d0 - probai(i))/slope
         endif
c   50% sup  
         if (probai(i).lt.0.75 .AND. probai(i+1).ge.0.75) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zpdzs(1)=x(i) + (0.75d0 - probai(i))/slope
         endif
c   68% inf  
         if (probai(i).lt.0.16 .AND. probai(i+1).ge.0.16) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zpdzi(2)=x(i) + (0.16d0 - probai(i))/slope
         endif
c   68% sup
         if (probai(i).lt.0.84 .AND. probai(i+1).ge.0.84) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zpdzs(2)=x(i) + (0.84d0 - probai(i))/slope
         endif
c   90% inf  
         if (probai(i).lt.0.05 .AND. probai(i+1).ge.0.05) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zpdzi(3)=x(i) + (0.05d0 - probai(i))/slope
         endif
c   90% sup  
         if (probai(i).lt.0.95 .AND. probai(i+1).ge.0.95) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zpdzs(3)=x(i) + (0.95d0 - probai(i))/slope
         endif
c   99% inf  
         if (probai(i).lt.0.005 .AND. probai(i+1).ge.0.005) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zpdzi(4)=x(i) + (0.005d0 - probai(i))/slope
         endif
c   99% sup  
         if (probai(i).lt.0.995 .AND. probai(i+1).ge.0.995) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zpdzs(4)=x(i) + (0.995d0 - probai(i))/slope
         endif
      enddo
c   
      return
      END
c
