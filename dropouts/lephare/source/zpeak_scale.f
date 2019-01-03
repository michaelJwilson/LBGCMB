c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Search for high peaks in The ML function vs z
c     And sort them from highest to smallest peaks in ML
      SUBROUTINE ZPEAK_SCALE(chio,y,nmax,zstep,scale,thres,zsec,dnz)
c     nmax = number of z 
c     dnz  = number of peaks in the Chi2 function
c     zsec(1-->dnz) number associate to the z 
      implicit none
      integer*4 chisize,npf
      INCLUDE 'dim_zchi.decl'
      parameter (npf=100)    
      REAL*8    chio(2,chisize),x(chisize),y(chisize)
      REAL*8    zpeak(0:chisize),max(0:chisize)
      REAL*8    scale,zsels(npf),mlsels(npf)
      REAL*8     compac(2,0:chisize),zstep,maxi
      REAL*8    thres,zsel(npf),mlsel(npf),max1,zpeak1
      INTEGER*4 I,j,k,nwin,nmax,indx(npf),id1
      INTEGER*4 dnz,step,kmax,zsec(npf)
      INTEGER*4 id(0:chisize),idcomp(0:chisize),idsel(npf)
c
c Initialisation
c      scale=0.2   ! scale in Dz for the window with search of maxima 
c      thres=0.1   ! level uses to considerate a peak detection  
c
      zpeak1= 0
      max1  = 0
      id1   = 0
      dnz   = 0
      do i = 0,chisize
        zpeak(i)   = 0
        max(i)     = 0
        compac(1,i)= 0
        compac(2,i)= 0
        id(i)      = 0
        idcomp(i)  = 0
      enddo
      do i = 1, npf
        zsels(i) = 0
        mlsels(i)= 0
        zsel(i)  = 0
        mlsel(i) = 0
        idsel(i) = 0
        zsec(i)  = 1
        indx(i)  = 0
      enddo
c               
      do i = 1, nmax
         x(i) = chio(1,i)
      enddo
c
      nwin=IDINT(DINT(scale/zstep))
      if (nwin.ge.IDINT(DBLE(nmax)/2.)) nwin=IDINT(DBLE(nmax)/2.-1.)
      if (nwin.eq.0) nwin=1
      step = 1 
c
c Search for all the maxima and initialize the extrema       
      k        = 0 
      zpeak(k) = -1
      max(k)   = 0.     
      do i=nwin,nmax-nwin,step
         maxi=0
         k=k+1
c inside the window 2*nwin +1 
         do j = i-nwin+1,i+nwin
            if (y(j) .gt. maxi) then
              zpeak1= x(j)
              max1  = y(j)
              id1   = j
              maxi  = y(j)
            endif
         enddo
         if (maxi.gt.0) then
            zpeak(k)= zpeak1
            max(k)  = max1
            id(k)   = id1
         else 
            zpeak(k)= x(i)
            max(k)  = maxi
            id(k)   = i
         endif            
      enddo
      kmax = k
      zpeak(kmax+1) = x(nmax)+1
      max(kmax+1)   = 0.
c Reducing to the k different values 
      k=1
      compac(1,k) = max(0)
      compac(1,k) = zpeak(0)
      do i=0,kmax
        if (max(i).ne.max(i+1)) then
          k = k + 1
          compac(1,k)= max(i) 
          compac(2,k)= zpeak(i)
          idcomp(k)  = id(i)
        endif
      enddo
      k=k+1
      compac(1,k) = max(kmax+1) 
      compac(2,k) = zpeak(kmax+1)
c Search inside these k values the maxima          
      dnz = 0
      do i = 2,k-1   
         if (compac(1,i).gt.thres) then
           if (compac(1,i).gt.compac(1,i-1).and.
     >       compac(1,i).gt.compac(1,i+1)) then
             dnz        = dnz+1
             zsel(dnz)  = compac(2,i)
             mlsel(dnz) = compac(1,i)
             idsel(dnz) = idcomp(i)
           endif
         endif
      enddo        
c  sort the z with ML values 
      if (dnz.gt.1) then
         call indexx(dnz,mlsel,indx)
         do k = 1,dnz
          zsels(dnz-k+1)  = zsel(indx(k))
          mlsels(dnz-k+1) = mlsel(indx(k))
          zsec(dnz-k+1)   = idsel(indx(k)) 
         enddo
      else
          zsels(1)  = zsel(1)
          mlsels(1) = mlsel(1)      
          zsec(1)   = idsel(1)
      endif
c
      return 
      end
