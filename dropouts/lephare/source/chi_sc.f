c       
c       program    CHI_SC
c
	subroutine CHI_SC(imodmax,fluxmod,distmod,tuniv,
     >  paramod,parasc_lim,iscmax,
     >  zobs,ab,sab,nbf,imagm,bused,parasc) 
c              
c   * Reads the S Charlot MODELs at various Z for fixed bands
c     as defined in SC models (BVRIJK)VVDS + (ugriz)CFHTLS + (IRAC1234)    
c     build the list of models Mag according to filters :
c     FILTER_FILE : /data/arnouts/zpwork/filt/vvds_cfhtls_swire.filt
c     U B V R I up gp rp ip zp J Ks IRAC_1 IRAC_2 IRAC_3 IRAC_4 24mic
c      Output format 
c      Model Z  Age Age_wr TauV mu  M*  M*/(1-R) SFR(1e8)
c      U B V R I u g r i z J K IRAC1 2 3 4 24mic 
        implicit none 
c
c  library 
      integer*4   sclib,iscmax,xsize,nmod
      PARAMETER   (sclib=90000)
      PARAMETER   (xsize=250)
      real*8    paramod(15,10,sclib),parasc_lim(2,10)     
c      real*8    zmod(15,sclib)
      real*8    tuniv,distmod
      real*8    fluxmod(15,17,sclib),fluxm(17)
      integer*4 imodmax(15),bdmax
c  observation 
      integer*4 imagm,nbf,bused(nbf)
      real*8    zobs,ab(nbf),sab(nbf)
c  outputs 
      real*8    parasc(10,3) 
c  internal
      integer*4 i,j,k,pass
      integer*4 nf
c      real*8    zage,dage
      real*8    dm,avmago,avmagt,chi2
      real*8    xpar(10,xsize),ypar(10,xsize),xpara
      real*8    xstep,xmin,xmax
      real*8    x(xsize),y(xsize),area
      real*8    modmed,modinf(10),modsup(10)


c  step size for the physical parameters  
      xstep=0.03       ! change also xsize
      bdmax=13
c
c      open(91,file='test.chi')
c      open(92,file='test2.chi')
c  initialisation of outputs
      do k = 1, 10
         do i = 1,3
            parasc(k,i) = -999
         enddo
      enddo   
c  build the chi2 and para-axis
      do k = 1,7 
         do i = 1, xsize
           if (k.lt.5)  xpar(k,i) = parasc_lim(1,k) + (i-1.5)*xstep
           if (k.eq.5 .or. k.eq.6 ) xpar(k,i) =  5. + (i-1)*xstep
           if (k.eq.7)              xpar(k,i) = -4. + (i-1)*xstep
           ypar(k,i) = 0.d0
         enddo
      enddo
c  get the closest z value from SC_LIBRARY    
c   for the moment from 0.1 to 1.5        
c      zage = 0.  
c      dage = 999.
c      do i = 1,15
c         if ( dabs(zobs-0.1*i).lt.0.1 .and.
c     >        dabs(zobs-0.1*i).lt.dage) then
c           zage = i*0.1
c           dage = dabs(zobs-zage)
c         endif
c      enddo
       nmod = NINT(10.*zobs) 
       if (nmod.lt.1)  nmod=1
       if (nmod.gt.15) nmod=15
     
c      write(*,*) ' '
c      write(*,'(A,1x,f8.3,1x,i4,2(1x,e12.6))') 'chi2 : ',
c     >     zobs,nmod,distmod,tuniv
c      write(*,*) ' '
c
c  START LOOP ON SC_LIBRARY 
      pass=0
      do i = 1,imodmax(nmod)
c  check if age older than age of universe 
         if (zobs.gt.1.5.and.dlog10(tuniv).lt.paramod(nmod,1,i)) goto 2 
c  measuring scaling factor 
           avmago = 0. 
           avmagt = 0. 
           nf = 0 
c           do k = 2,14  ! B->4.5mic
           do k = 2,bdmax  ! B->3.6mic
              fluxm(k) = fluxmod(nmod,k,i)*10**(-0.4*distmod)
              if (ab(k).gt.0.and.sab(k).gt.0.and.
     >            bused(k).eq.1) then
                nf = nf + 1
       avmago = avmago + ab(k)*fluxm(k)/sab(k)**2
       avmagt = avmagt + fluxm(k)*fluxm(k)/sab(k)**2
              endif
           enddo
           if (nf.eq.0 .and. i.eq.1) then
             goto 1          ! exit 
           elseif (nf.eq.1 .and. i.eq.1 ) then 
             goto 1           ! exit
           else   
             if (avmagt.lt.1.d-60) goto 2
             dm = avmago/avmagt
           endif      
c  measuring chi2
           chi2=0.d0
c           do k = 2, 14
           do k = 2, bdmax
             if (bused(k).eq.1  .and. sab(k).gt.0) then
               chi2 = chi2 + ((ab(k)-dm*fluxm(k))/sab(k))**2
             endif
           enddo           
           if (chi2.lt.1.e9) chi2 = dexp(-0.5*chi2)
c  store chi2 vs para value
           do k = 1, 7 
              xpara = paramod(nmod,k,i) 
              if (k.ge.5)  xpara=paramod(nmod,k,i) + dlog10(dm)
c           if (k.eq.5 .and. MOD(i,10000).eq.0) write(*,*) dm,xpara
              do j = 1, xsize
                 if (DABS(xpara-xpar(k,j)).le.(xstep/2.)) then         
                    ypar(k,j)= ypar(k,j)+chi2
                 endif 
              enddo
           enddo
c
 2    continue
c  close big loop   
      enddo    ! endif sc_lib big loop              


c
c  estimates of the baysian parameters
      do k = 1,7
        xmin = xpar(k,1)
        xmax = xpar(k,xsize)
c        write(*,*) xmin,xmax
        do i = 1,xsize
           x(i) = xpar(k,i)
           y(i) = ypar(k,i)
c           write(91,*) k,x(i),y(i)
c           if (y(i).ne.0) write(*,*) i,x(i),y(i)
        enddo
c        call TRAPZD(x,y,xsize,xmin,xmax,area)
        if (k.eq.1) call TRAPZD(x,y,xsize,xmin,xmax,area)
c        if (k.eq.1 .or.k.eq.7) call TRAPZD(x,y,xsize,xmin,xmax,area)
c        write(*,*) k,area
        if (area.gt.0) then 
           call PROBAZBAY(x,y,xsize,area,modmed,modinf,modsup)
           parasc(k,1) = modmed
           parasc(k,2) = modinf(2)
           parasc(k,3) = modsup(2)
        endif     
c        write(*,*) k,(parasc(k,i),i=1,3)
      enddo

c      close(91)


 1    return


      end


c
