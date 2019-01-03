c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Return zmin,zmax included inside 50%,68%,90%,95% of the area 
c     The Pdz curve peaks at zbest Pdz=1
c
c      PROGRAM PROBAZ
      SUBROUTINE PROBAZ(x,y,bmax,area,zpdzi,zpdzs,pdzlevel)
      implicit none 
c      
      INTEGER*4 bmax,i,j,k,pass,step,taille,sens,l
      REAL*8    x(bmax),y(bmax),yn(bmax),a
      REAL*8    x2(10000),y2(10000),y3(10000)
      REAL*8    zli,zls,zi(bmax),zs(bmax)
      REAL*8    thres,dthres,area,sum,proba
      REAL*8    probai(1000),zinf(1000),zsup(1000)
      REAL*8    pdzlevel(10),zpdzi(10),zpdzs(10)

c
c   Initialisation
      pass  = 0
      zli   = 0.d0
      zls   = 0.d0
      proba = 0.d0
      do i = 1, bmax
         yn(i) = y(i)/area
         zi(i) = -1.d0
         zs(i) = -1.d0
      enddo      

c     Reinterpole la PDFz avec des pas 10 fois plus petits
      taille=(bmax-1)*10
      do i =1,taille
        x2(i)= x(1) + dble(i)*(x(bmax)-x(1))/dble(taille)
        l=int(dble(i)/10.)+1
        sens=idnint(dsign(1.d0,(x2(i)-x(l))))     
        a=dabs((x(l)-x2(i))/(x(l)- x(l+sens)))
        y2(i) = a* y(l+sens)+(1-a)*y(l)
        y3(i) = a*yn(l+sens)+(1-a)*yn(l)
      enddo
c
      dthres = 0.001d0      ! dthres > 0.01
      step = IDNINT(1.d0/dthres) 
      do i = 1, step
         probai(i)  = 0.d0
         zinf(i)   = -1.d0
         zsup(i)   = -1.d0
      enddo   
c
      do i = 1, 10
        pdzlevel(i) = 0.d0
        zpdzi(i) = -1.d0
        zpdzs(i) = -1.d0
      enddo  
c
c  thresholding the y-axis 
      do j = 1, step
        thres= 1.0 - DBLE(j)*dthres
        if (thres.gt.0) then 
c
        sum  = 0.   ! area below the curve at a fixed threshold between couples (zinf-zsup) 
        k    = 0    ! number of redshift interval used 
        pass = 0    ! compute area for each interval zinf-zsup
c  check Zinf at the limit
        if ( y2(1).ge.thres ) then
           k     = k + 1
           zi(k) = x2(1)
        endif
c  search for Zinf excluding boundaries 
        do i = 1, taille-1
           if (y2(i).lt.thres .AND. y2(i+1).gt.thres) then
              k = k + 1 
              zi(k) = x2(i) 
           endif
           if (y2(i).gt.thres .AND. y2(i+1).lt.thres) then 
              zs(k) = x2(i+1) 
              pass = 1
           endif
           if (i.eq.(taille-1) .AND. y2(i+1).ge.thres) then 
              zs(k) = x2(i+1) 
              pass = 1
           endif   
c  compute area between zinf-zsup 
           if (pass.eq.1) then
              zli = zi(k)
              zls = zs(k)
              call TRAPZD(x2,y3,taille,zli,zls,proba)
c              probai(j) = probai(j) + proba
              probai(j) = probai(j) + proba
              pass = 0
           endif     
c           write(*,*)"LA",i,thres,x2(i),y2(i),probai(j)          
        enddo   ! loop in z
c         write(*,*)j,thres,probai(j)

c
        if (probai(j) .gt. 0) then 
           zinf(j) = zi(1)            ! take the first Zinf 
           zsup(j) = zs(k)            ! take the last  Zsup
        endif   
c
        endif   ! threshold > 0 
c           
      enddo     ! loop on threshold ...
c
c  selection of the best Proba interval  
      do j = 1, step
         if (probai(j).le. 0.68 .AND. probai(j+1).gt.0.68) then
            zpdzi(1) = zinf(j)
            zpdzs(1) = zsup(j)
            pdzlevel(1) = (probai(j)+probai(j+1))/2.*100
         elseif  (probai(j).le. 0.90 .AND. probai(j+1).gt.0.90) then  
            zpdzi(2)    = zinf(j)
            zpdzs(2)    = zsup(j)
            pdzlevel(2) = (probai(j)+probai(j+1))/2.*100
         elseif  (probai(j).le. 0.99 .AND. probai(j+1).gt.0.99) then  
            zpdzi(3)    = zinf(j)
            zpdzs(3)    = zsup(j)
            pdzlevel(3) = (probai(j)+probai(j+1))/2.*100
         endif
      enddo



c   
      return
      END
c
