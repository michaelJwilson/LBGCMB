c      last modif 06/01/2000
c      Author : S. ARNOUTS 
c      description : It allows to resample two initial functions
c      (xi1,yi1) length=imax1,(xi2,yi2) length=imax2
c      with the same abscisse ls (length=smax)
c      in the common range of (xi1,xi2) and to interpolate the 
c      ordinates y1 and y2 according to the new abscisse (ls)  
c      Input : (xi1,yi1),imax1  and (xi2,yi2), imax2
c      Output : (ls,y1) and (ls,y2) , smax 
c
      subroutine sampling(x1i,y1i,imax1,x2i,y2i,imax2,ls,y1,y2,smax)
      implicit none 
c  input 
      integer*4  imax1,imax2
      real*8     x1i(imax1),y1i(imax1)
      real*8     x2i(imax2),y2i(imax2)
c  interm      
      integer*4 maxsize,wmax
      INCLUDE 'dim_wave.decl'
      parameter (maxsize=110000)
c
c      integer*4 max1,k,i,j,iwksp(maxsize),ktest
c      real*8    xmin,xmax,undefval
c      real*8    y1s(maxsize),y2s(maxsize) 
c      real*8    xt(maxsize),yt1(maxsize),yt2(maxsize) 
c      real*8    sav1(maxsize),sav2(maxsize),sav3(maxsize)
c      real*8    xt10,yt10,xt1n,yt1n,xt20,yt20,xt2n,yt2n
c output
c      real*8     ls(maxsize),y1(maxsize),y2(maxsize)
c      integer*4  smax
c
      integer*4 max1,k,i,j,iwksp(wmax)
      real*8    xmin,xmax,undefval
      real*8    y1s(wmax),y2s(wmax) 
      real*8    xt(wmax),yt1(wmax),yt2(wmax) 
      real*8    sav1(wmax),sav2(wmax),sav3(wmax)
      real*8    xt10,yt10,xt1n,yt1n,xt20,yt20,xt2n,yt2n
c output
      real*8     ls(wmax),y1(wmax),y2(wmax)
      integer*4  smax
c
      parameter (undefval=-999.999)
c
      do i = 1, wmax
         ls(i) =0 
         y1(i) =0 
         y2(i) =0 
         xt(i) =0 
         yt1(i) =0 
         yt2(i) =0 
      enddo
c      write(*,*) imax1,imax2 
c  Search the common range in abscisse
      xmin=0.d0
      xmax=0.d0
      if (x1i(1).ge.x2i(1)) xmin=x1i(1)
      if (x1i(1).lt.x2i(1)) xmin=x2i(1) 
      if (x1i(imax1).ge.x2i(imax2)) xmax=x2i(imax2)
      if (x1i(imax1).lt.x2i(imax2)) xmax=x1i(imax1)
c      write(*,*) xmin,xmax 
c  Redefined a unique vector xt inside range
      xt10=undefval
      yt10=undefval
      xt1n=undefval
      yt1n=undefval
      xt20=undefval
      yt20=undefval
      xt2n=undefval
      yt2n=undefval
      k = 0
      do i = 1,imax1
         if (i.lt.imax1.and.((x1i(i)-xmin)*(x1i(i+1)-xmin)).lt.0.d0)then
           xt10=x1i(i)
           yt10=y1i(i)
         endif
         if (i.gt.1.and.((x1i(i-1)-xmax)*(x1i(i)-xmax)).lt.0.d0) then
           xt1n=x1i(i)
           yt1n=y1i(i)
         endif      
         if (x1i(i).ge.xmin .and. x1i(i).le.xmax) then
            k = k + 1
            xt(k) = x1i(i)
            yt1(k)= y1i(i)
            yt2(k) = undefval 
         endif
      enddo
      max1=k
      do i = 1,imax2 
         if (i.lt.imax2.and.((x2i(i)-xmin)*(x2i(i+1)-xmin)).lt.0.d0)then
           xt20=x2i(i)
           yt20=y2i(i)
         endif
         if (i.gt.1.and.((x2i(i-1)-xmax)*(x2i(i)-xmax)).lt.0.d0) then
           xt2n=x2i(i)
           yt2n=y2i(i)
         endif    
         if (x2i(i).ge.xmin .and. x2i(i).le.xmax) then
c  check if xt exists already ?
            do j = 1, max1
              if (DABS(x2i(i)-xt(j)).le.1.d-9) then
                yt2(j)=y2i(i)  
                goto 19
              endif
            enddo
            k = k + 1
            xt(k)= x2i(i)
            yt2(k)=y2i(i)
            yt1(k) = undefval
 19         continue
         endif
      enddo      
      smax=k
c      ktest=0
c      do k = 1,smax
c         xtt(k)=xt(k)
c      enddo
c      write(*,*) 'ok ',smax
c  sort with respect to the abscisses
      call indexx(smax,xt,iwksp) 
c      write(*,*) 'ok ',smax,xt(iwksp(1)),xt(iwksp(smax))
      do i = 1,smax
         sav1(i) = xt(i)
         sav2(i) = yt1(i)
         sav3(i) = yt2(i)
      enddo
      do i = 1, smax
         ls(i) = sav1(iwksp(i))
         y1s(i) = sav2(iwksp(i))
         y2s(i) = sav3(iwksp(i))
      enddo
c  Interpolation for missing value 
      call interpol(ls,y1s,y2s,smax,xt10,yt10,xt1n,yt1n,
     $xt20,yt20,xt2n,yt2n,y1,y2)
c      do i = 1,smax 
c        if (y1s(i).eq.undefval) y1(i)=0.d0 
c        if (y2s(i).eq.undefval) y2(i)=0.d0 
c      enddo
c    
c      write(*,*) imax1,imax2,smax
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE interpol(ls,y1s,y2s,smax,x10,y10,x1n,y1n,
     $x20,y20,x2n,y2n,y1,y2)
      implicit none 
      integer*4 maxsize,wmax
       INCLUDE 'dim_wave.decl'    
      parameter (maxsize=110000)
c      
c      integer*4 smax,i,k,kinf,ksup,smax2
c      real*8    undefval,ls(maxsize),y1s(maxsize),y2s(maxsize)
c      real*8    y1(maxsize),y2(maxsize),y1c(0:maxsize),y2c(0:maxsize)
c      real*8    x10,y10,x1n,y1n,x20,y20,x2n,y2n,lsc(0:maxsize)
c      parameter (undefval=-999.999)
c
      integer*4 smax,i,k,kinf,ksup,smax2
      real*8    undefval,ls(wmax),y1s(wmax),y2s(wmax)
      real*8    y1(wmax),y2(wmax),y1c(0:wmax),y2c(0:wmax)
      real*8    x10,y10,x1n,y1n,x20,y20,x2n,y2n,lsc(0:wmax)
      parameter (undefval=-999.999)
c
c      write(*,*) x10,y10,x1n,y1n,x20,y20,x2n,y2n
      smax2 = smax + 1
      if (x10.ne.undefval) then
        y1c(0)=y10
        lsc(0)=x10
      endif
      if (x1n.ne.undefval) then
        y1c(smax2)=y1n
        lsc(smax2)=x1n
      endif
      if (x20.ne.undefval) then
        y2c(0)=y20
        lsc(0)=x20
      endif
      if (x2n.ne.undefval) then
        y2c(smax2)=y2n
        lsc(smax2)=x2n
      endif
      do i = 1, smax
        y1c(i)=y1s(i)
        y2c(i)=y2s(i) 
        lsc(i)=ls(i)
      enddo
      do i = 1, smax
        y1(i)=y1c(i) 
        y2(i)=y2c(i) 
        if (y1(i).eq.undefval) then
          k = i
          do  while (y1c(k).eq.undefval) 
             k = k - 1
          enddo
          kinf=k
          k = i
          do  while (y1c(k).eq.undefval) 
             k = k + 1
          enddo
          ksup=k
          y1(i)=y1c(kinf)+(lsc(i)-lsc(kinf))*
     $         (y1c(ksup)-y1c(kinf))/(lsc(ksup)-lsc(kinf))  
        endif
        if (y2(i).eq.undefval) then
          k = i
          do  while (y2c(k).eq.undefval) 
             k = k - 1
          enddo
          kinf=k
          k = i
          do  while (y2c(k).eq.undefval) 
             k = k + 1
          enddo
          ksup=k
          y2(i)=y2c(kinf)+(lsc(i)-lsc(kinf))*
     $         (y2c(ksup)-y2c(kinf))/(lsc(ksup)-lsc(kinf)) 
        endif
      enddo
 
      return
      end
c



