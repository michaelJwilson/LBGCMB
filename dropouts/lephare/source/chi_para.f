      subroutine chi_para(zchipara,liblength,zinf,zsup,zstep,
     >       parainf,parasup,paramed)  
c 
      implicit none 
c
      integer*4   inlib
      INCLUDE    'dim_lib.decl'
c
      integer*4  liblength,i,j,k,lmax(6),lirmaxs,nchi
      integer*4  valtest,test
      real*8     zchipara(12,inlib),parainf(6),parasup(6),paramed(6) 
      real*8     zinf,zsup,zstep
      real*8      bstart,bend,barea,dbin(6),sbin(6),ebin(6)
      real*8     physp(6,300),chip(6,300)
      real*8     lumir(300),chilir(300)
      real*8     zmed,zpdzi,zpdzs
c
c     zchipara :
c        1  2   3    4    5     6    7   8   9      10   11  12   
c        z,age,ext-l,ebv,ldust, luv, lr, lk, ldustp, mo, sfr  chi2
c        N   Y   N    N   Y     N    N    N    Y     Y    Y   N       
c        <--  mag_gal file --> <----    phys file        --->

c      only use : age ldust ldustp m sfr + ssfr

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  number of resolution elements per phys. parameters       
      do i = 1, 6
        dbin(i)=0
        ebin(i)=0
        sbin(i)=0
        lmax(i)=101         
      enddo
c       
c  log Age  [yr]  ! 7 -> 10.5
      dbin(1) = 0.05
      sbin(1) = 7.0
      lmax(1) = 71
c  Ldust [L0]  ! from mag_gal (EB-V, Extlaw)  6 to 14
      dbin(2) = 0.05                 
      sbin(2) = 6.0
      lmax(2) = 161
c  Ldust [L0]  ! from  sedtolib (*.phys)  6 to 14 
      dbin(3) = 0.05                 
      sbin(3) = 6.0
      lmax(3) = 161
c  log Stellar Mass  [Mo] : 7 to  14
      dbin(4) = 0.05
      sbin(4) = 7.
      lmax(4) = 151
c log SFR 1.e8/1.e9  [Mo/yr] : -6 to 6
      dbin(5) = 0.1
      sbin(5) = -6.0
      lmax(5) = 121
c log SFR/Mass [SFR/Mo] : -17.5 to -6.5 
      dbin(6) =  0.1
      sbin(6) = -17.5
      lmax(6) =  111 
cccccccc
      do j = 1, 6
        do k = 1, lmax(j) 
           physp(j,k)= sbin(j)+dbin(j)*(k-1)  
           chip(j,k) = 0.d0
        enddo
        ebin(j)=sbin(j)+dbin(j)*(lmax(j)-1)  
c        write(*,*) physp(j,1) , physp(j,lmax(j))
      enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   loop on library 
      do i = 1, liblength 
        if (zchipara(1,i).ge.(zinf-zstep)
     >    .and. zchipara(1,i).le.(zsup+zstep)
     >    .and. zchipara(12,i).le.1.e8
     >    .and. zchipara(12,i).gt.0) then
c          write(*,*) zchipara(2,i),zchipara(5,i),zchipara(9,i),
c     >     zchipara(10,i),zchipara(11,i),zchipara(11,i)-zchipara(10,i)

c log  age
          do nchi = 1,lmax(1)                    
            if (SNGL(DABS(dlog10(zchipara(2,i))-physp(1,nchi)))
     >                                 .le.SNGL(dbin(1)/2.d0))then
             chip(1,nchi)=chip(1,nchi)+dexp(-0.5*zchipara(12,i))
             goto 1
c        write(*,*) dlog10(zchipara(2,i)),physp(1,nchi),
c     >            dbin(1),chip(1,nchi)
            endif
          enddo
c log ldust
 1        do nchi = 1,lmax(2)                    
           if (SNGL(DABS( zchipara(5,i)-physp(2,nchi) ))
     >                            .le.SNGL(dbin(2)/2.d0) ) then
            chip(2,nchi)=chip(2,nchi)+dexp(-0.5*zchipara(12,i))
             goto 2
           endif
          enddo
c log ldust
 2        do nchi = 1,lmax(3)                     
           if (SNGL(DABS( zchipara(9,i)-physp(3,nchi) ))
     >                             .le.SNGL(dbin(3)/2.d0) ) then
            chip(3,nchi)=chip(3,nchi)+dexp(-0.5*zchipara(12,i))
             goto 3
           endif  
          enddo
c log mass
 3        do nchi = 1,lmax(4)                    
           if (SNGL(DABS( zchipara(10,i)-physp(4,nchi) ))
     >                               .le.SNGL(dbin(4)/2.d0) )then
            chip(4,nchi)=chip(4,nchi)+dexp(-0.5*zchipara(12,i))
             goto 4
           endif
          enddo
c log  sfr
 4        do nchi = 1,lmax(5)                    
          if (SNGL(DABS( zchipara(11,i)-physp(5,nchi) ))
     >                               .le.SNGL(dbin(5)/2.d0) )then
             chip(5,nchi)=chip(5,nchi)+dexp(-0.5*zchipara(12,i))
              goto 5
           endif
          enddo
c  log SSFR1e9
 5        do nchi = 1,lmax(6)                    
            if(SNGL(DABS( zchipara(11,i)-zchipara(10,i)
     >           -physp(6,nchi) )).le.SNGL(dbin(6)/2.d0) )then
           chip(6,nchi)=chip(6,nchi)+dexp(-0.5*zchipara(12,i))
               goto 6
            endif
          enddo
 6        continue
        endif
      enddo 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   measuring the median and quartiles 
c      write(UO,*) ' analysis 0'
      valtest=-99    
      do k = 1,6
        paramed(k) = -99.
        parainf(k) = -99.   
        parasup(k) = -99.  
        test=0
        do j = 1,lmax(k)
           if (chip(k,j).lt.1.d-60) chip(k,j)=0.d0
           lumir(j) =physp(k,j)
           chilir(j)=chip(k,j)
           if (chilir(j).gt.0) test = 1
        enddo
c        write(*,*) 'test',test
        if (test.gt.0) then 
            bstart= sbin(k)
            bend  = ebin(k)
            lirmaxs=lmax(k)
            call TRAPZD(lumir,chilir,lirmaxs,bstart,bend,barea)
c        write(*,*) test,barea
            if (barea.gt.0) then 
           call PROB_MED2(lumir,chilir,lirmaxs,barea,zmed,zpdzi,zpdzs)
               paramed(k) = zmed
               parainf(k) = zpdzi   
               parasup(k) = zpdzs     
            endif
            if (k .eq.1) then 
               paramed(1) = 10**(paramed(1))
               parainf(1) = 10**(parainf(1))
               parasup(1) = 10**(parasup(1))
            endif
c            write(*,*) parainf(k),paramed(k),parasup(k)
        endif
      enddo 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Return zmin,zmax included inside 50%,68%,90%,95% of the area 
c
c      PROGRAM PROBAZBAY
      SUBROUTINE PROB_MED2(x,y,bmax,area,zmed,zpdzi,zpdzs)
      implicit none 
c      
      INTEGER*4 bmax,i
      REAL*8    x(300),y(300),yn(300),zmed,slope
      REAL*8    zli,zls,area,proba,probai(300),zpdzi,zpdzs
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
      zpdzi = -99.
      zpdzs = -99.
cccc Method with Area  met=0
        probai(1) = 0.d0
        do i = 1, bmax-1
           zls=x(i+1)
           call TRAPZD(x,yn,bmax,zli,zls,proba)
           probai(i+1) = proba
        enddo
        do i = 1, bmax-1
c   Val median 
          if (probai(i).lt.0.5d0 .AND. probai(i+1).ge.0.5d0) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zmed=x(i) + (0.5d0 - probai(i))/slope
  
c             if (DABS(probai(i+1)-0.5d0).lt.0.00001d0) then
c               zmed = x(i+1) 
c             else
c               zmed = (x(i+1)+x(i))/2.d0 
c             endif
          endif         

c   68% inf  
          if (probai(i).lt.0.16d0 .AND. probai(i+1).ge.0.16d0) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zpdzi = x(i) + (0.16d0 - probai(i))/slope

c             if (DABS(probai(i+1)-0.16d0).lt.0.00001d0) then
c                zpdzi = x(i+1)
c             else
c                zpdzi = (x(i+1)+x(i))/2.d0 
c             endif
          endif

c   68% sup
          if (probai(i).lt.0.84d0 .AND. probai(i+1).ge.0.84d0) then
             slope= ( probai(i+1)-probai(i) )/( x(i+1)-x(i) )
             zpdzs=x(i) + (0.84d0 - probai(i))/slope

c             if (DABS(probai(i+1)-0.84d0).lt.0.00001d0) then
c                zpdzs = x(i+1)
c             else 
c                zpdzs = (x(i+1)+x(i))/2.d0 
c             endif

          endif
        enddo

      return
      END
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
