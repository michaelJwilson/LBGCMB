      subroutine chi_libphys(libphys,libppara,
     >      zrecp,zrecpi,zrecps,nzrecp,isubmax,
     >      zphot,ab,sab,bused,buscal,lirmed,
     >      ppbest,chipbest,reclphysb,dmpbest,ppmed,ppinf,ppsup,nppara,
     >      fluxphys,kcorphys,magphys0)  

      implicit none 
      integer*4 chisize,nbf,maxsize
      parameter (maxsize=110000)  
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_zchi.decl'
      INCLUDE 'out_unit.decl' 
      integer*4     imag,test,vali,nchi,npsum,iter
      integer*4     i,j,k,irec,iz,lirmax(50),lirmaxs
      integer*4     reclmax,modlib,reclib
      character*4096 param,paravc(500),file,zpdir,zpwork
      real*8         physpmin(50,500),chiphys(50,500)
      real*8         lumir(500),chilir(500),chipbest,lirmed
c      real*8         phys_cont,phys_scale   
c      integer*4     nbuphys,nbsphys,busphys(nbf),bscphys(nbf)  
      integer*4      nf_phys,reclphysb,isub,isubmax
      real*8         dm,avmagt,avmago,chi2_phys,chirmin,chiorig
      real*8        dz,val,val2,h0,om0,l0,funz,zlib,dmcor
      real*8        bstart,bend,dbin(50),sbin(50),ebin(50)
      real*8        barea,zmed,zpdzi,zpdzs,valtest,ztemp,mod_dist
      real*8        facerr,sabo(nbf)
c
c input 
      character*4096 libphys
      real*8        libppara(50,maxsize),zphot,tuniv,timy
      integer*4     bused(nbf),buscal(nbf)   
      integer*4     nppara(350)
      real*8        ab(nbf),sab(nbf),maglib(nbf),maglibf(nbf)
      real*8        zrecp(500),kcorlib(nbf)
      integer*4     zrecpi(500),zrecps(500),nzrecp
c output 
      real*8        ppbest(50),ppmed(50),ppinf(50),ppsup(50)
      real*8        fluxphys(nbf),dmpbest
      real*8        kcorphys(nbf),magphys0(nbf)
c
      real*8      Lsol
      parameter  (Lsol=3.826e33)
      external      funz,timy
ccccccccccccccccccccc
      if(UO.ne.6)
     >   open(UO,file='/tmp/screenChiPhyslib.dat',status='unknown') 
      call getenv('LEPHAREDIR',zpdir)
      if (lnblnk(zpdir) .eq. 0) then
          write(UO,*) 'WARNING :  variable LEPHAREDIR not defined'
          stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      if (lnblnk(zpwork) .eq. 0) then
          write(UO,*) 'WARNING :  variable LEPHAREWORK not defined'
          stop
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   extracting  the record for specific z
c       
      file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >            // libphys(1:lnblnk(libphys)) // '.doc'
      param='RECORD_LENGTH'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') reclmax
      param='Z_STEP'
      call read_para2(param,file,paravc,test)
      read(paravc(1),'(f8.3)') dz
      param='COSMOLOGY'
      call read_para2(param,file,paravc,test)
      read(paravc(1),'(f8.3)') h0
      read(paravc(2),'(f8.3)') om0
      read(paravc(3),'(f8.3)') l0
c
      iz=0
      dmcor=1.d0
      if( SNGL(DABS(zphot-zrecp(1))).le.
     >        SNGL((zrecp(2)-zrecp(1))/2.d0) ) then 
          iz=1
          ztemp=zrecp(1) 
          tuniv  = timy(ztemp,h0,om0,l0)      ! 0 <-- zinfinity
          mod_dist=funz(ztemp,h0,om0,l0)
          dmcor=10**(0.4*(mod_dist-funz(zphot,h0,om0,l0)))
          goto 1
      endif
      do i = 2, nzrecp        
        if ( SNGL(DABS(zphot-zrecp(i))).le.
     >          SNGL((zrecp(i)-zrecp(i-1))/2.d0) ) then 
           iz=i
           ztemp = zrecp(i)-(zrecp(i)-zrecp(i-1))/2.d0
           tuniv = timy(ztemp,h0,om0,l0)      ! zlim <- zinfinity
           ztemp=zrecp(i)
           mod_dist=funz(ztemp,h0,om0,l0)
           dmcor=10**(0.4*(mod_dist-funz(zphot,h0,om0,l0)))
           goto 1
        endif 
      enddo
 1    if (iz.eq.0) goto 2
c      write(*,'(i8,1x,A)') iz,file(1:lnblnk(file))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  number of resolution elements per phys. parameters       
      iter=0
c      lirmax=101    
      npsum=0
      do i = 1, 50
        npsum=npsum+nppara(i)
        dbin(i)=0
        ebin(i)=0
        sbin(i)=0
        lirmax(i)=101         
      enddo
c       
      chiorig=0
      nf_phys  =-1
      dm       =-1
      chi2_phys=1.e9 
      chirmin=1.e9
      reclphysb=-1
c    parameters from sedtolib   
c          write(5,700) nrec,jstoch,tform 1,luv 2,lopt 3,lnir 4,ltir 5,mstr1 6 
c     >   (sfrav(i),i=3,4) 7,8 ,gamma 9,zmet 10 ,tauv0 11 ,mu 12,age_wm 13 ,age_wr 14
c     >    tlastburst 15,(fburst(i),i=3,4) 16 17 ,DBLE(nburst) 18 ,mstr0 19
c     >    lfuv 20 ,lnuv 21 ,loptr 22,lnirk 23,afuv 24 ,anuv 25 ,ag 26, ar 27 ,az 28 ak 29 ,ldusth 30, 
c     >    d4000 31

      if (npsum.gt.0) then 
        do j = 1,28
          if (j.eq.1 .OR. (j.ge.13.and.j.le.15)) then 
c                                                ! log Age  [yr]  ! log tlastburst [yr]
            dbin(j)   = 0.05                     !  5.5 to 10.5
            sbin(j)   = 5.5
            lirmax(j) = 105
         endif
          if (j.ge.2 .and. j.le.5) then          ! Log Luv , Lr , Lk , Ldust [Lo] :  6 to 14.5
            dbin(j)   = 0.05                 
            sbin(j)   = 6.0
            lirmax(j) = 171
          endif  
          if (j.eq.6) then                       ! log Stellar Mass  [Mo] :  6.5 to  14
            dbin(6)   = 0.03
            sbin(6)   = 6.5
            lirmax(6) = 251
          endif
          if (j.eq.7 .OR. j.eq.8) then           ! log SFR 1.e8/1.e9  [Mo/yr] : -12 to 8
            dbin(j)   = 0.05
            sbin(j)   = -12.0
            lirmax(j) = 401
          endif  
          if (j.eq.9) then                       ! log Age  [Gyr] : -4 to 1.2
            dbin(9)   = 0.05
            sbin(9)   = -4.0
            lirmax(9) = 105
          endif
          if (j.eq.10) then                      ! metallicity 0-2
            dbin(10)  = 0.05
            sbin(10)  = 0.0
            lirmax(10)= 41 
          endif
          if (j.eq.11) then                      ! tauv : 0 - 6  
            dbin(11)  = 0.05
            sbin(11)  = 0.0
            lirmax(11)= 121
          endif
          if (j.eq.12) then                      ! mu: 0 - 1. 
            dbin(12)  = 0.02
            sbin(12)  = 0.0
            lirmax(12)= 51 
          endif 
          if (j.eq.16 .OR. j.eq.17) then         ! frac burst over last 1.e8yr/1.e9yr : 0 to 1
            dbin(j)   = 0.01
            sbin(j)   = 0.0
            lirmax(j) = 101
          endif
          if (j.eq.18) then                      ! Nb of past bursts  : 0 to 20 
            dbin(18)  =  1.0
            sbin(18)  =  0.0
            lirmax(18)=  21
          endif
          if (j.eq.19 .OR. j.eq.20 ) then        ! log SFR/Mass 1.e8/1e9 [SFR/Mo] : -17 to -7 
            dbin(j)   =  0.05
            sbin(j)   = -17.0
            lirmax(j) =  201 
          endif         
          if (j.ge.21 .and.j.le.26) then         ! Attenuation Afuv,nuv,g,r,z,k {mag] : 0 to 10
            dbin(j)   = 0.04
            sbin(j)   = 0.0
            lirmax(j) = 250  
          endif
          if (j.eq.27) then                      !  Ldust HOT  (from Absorption with lbda<4000A)
            dbin(j)   = 0.05                 
            sbin(j)   = 6.0
            lirmax(j) = 171
          endif 
          if (j.eq.28) then                      !  D4000 : 0 to 3 
            dbin(j)   = 0.02                 
            sbin(j)   = 0.7
            lirmax(j) = 100
          endif 

cccccccc
           do k = 1, lirmax(j) 
              physpmin(j,k)= sbin(j)+dbin(j)*(k-1)  
              chiphys(j,k) = 0.d0
           enddo
           ebin(j)=sbin(j)+dbin(j)*(lirmax(j)-1)  
        enddo
      endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   loop on library 
 3    continue
      do k = 1,nbf
         fluxphys(k) = 0.d0
         kcorphys(k) = 0.d0 
         magphys0(k) = 0.d0 
      enddo
      isub=0
      file = zpwork(1:lnblnk(zpwork)) // '/lib_mag/'
     >       // libphys(1:lnblnk(libphys)) // '.bin'
c
      open(29,file=file(1:lnblnk(file)),status='unknown',
     >       access='direct',recl=reclmax)
      do irec = zrecpi(iz),zrecps(iz)
         read(29,rec=irec) modlib,vali,val,zlib,
     >        val2,reclib,imag,(maglib(k),k=1,imag),
     >       (kcorlib(k),k=1,imag)     
         if (SNGL(DABS(zlib-zphot)).gt.SNGL(dz/2.d0)) 
     >       write(UO,*) 'wrong redshift',zlib,zphot
c         write(*,*) irec,modlib,reclib,libppara(10,reclib)
         isub=isub+1
         if (isub.gt.isubmax)  goto 4 
c   Check that Age-Massweighted is smaller than Age of Universe
         if (libppara(13,reclib).gt.tuniv) goto 6 
         dm  = 0.
         avmagt= 0. 
         avmago= 0. 
         nf_phys=0 
         do k = 1,imag
            maglibf(k) = 10**(-0.4*(maglib(k)+48.59))
            if (ab(k).gt.0.and.sab(k).gt.0.and.maglibf(k).gt.0
     >         .and. bused(k).eq.1 .and.buscal(k).eq.1 ) then
             avmago=avmago+ab(k)*maglibf(k)*dmcor/sab(k)**2
             avmagt=avmagt+maglibf(k)*maglibf(k)*dmcor**2/sab(k)**2
                nf_phys = nf_phys + 1                     
            endif
         enddo
         if (nf_phys.ge.1) then
             dm = avmago/avmagt
         else
             dm = 1.d0
         endif
         chi2_phys = 0.d0
c  reject models with Ldust that disagree with L_FIR from FIR analysis (if measured)  
c         if (lirmed.gt.3 .and. lirmed.lt.17) then 
c          if (ABS(lirmed-libppara(5,reclib)-dlog10(dm)).gt.1) goto 6 
c          if ((lirmed-libppara(5,reclib)-dlog10(dm)).gt. 5.0) goto 6 
c          if ((lirmed-libppara(5,reclib)-dlog10(dm)).lt.-5.0) goto 6 
c         endif 
         nf_phys=0
         do k = 1, imag
            if (sab(k).gt.0.and.bused(k).eq.1) then
          chi2_phys=chi2_phys+((ab(k)-dm*maglibf(k)*dmcor)/sab(k))**2
               nf_phys=nf_phys+1
            endif 
         enddo
cc   RETURN  if less than 3 filters used
         if (nf_phys .le. 3) then 
            close(29) 
            goto 2
         endif
ccc   best chi2
         if (chi2_phys.lt.chirmin ) then
            chirmin   = chi2_phys
            dmpbest   = dm
            chipbest  = chi2_phys
            reclphysb = reclib 
            ppbest(1) = 0.d0
c            write(*,*) irec,modlib,reclib 
            do k = 1, imag
               fluxphys(k)=dm*maglibf(k)*dmcor
               kcorphys(k) = kcorlib(k)
            enddo
            if (libppara(1,reclib).gt.0) 
     >      ppbest(1) = dlog10(libppara(1,reclib))                  ! Lg Age
c            ppbest(2) = libppara(2,reclib)+dlog10(3.d18*400/2300**2)
c            ppbest(2) = ppbest(2)+dlog10(dm)-dlog10(Lsol)           ! log Luv
c            ppbest(3) = libppara(3,reclib)+dlog10(3.d18*1000/6000**2)
c            ppbest(3) = ppbest(3)+dlog10(dm)-dlog10(Lsol)           ! log Lr
c            ppbest(4) = libppara(4,reclib)+dlog10(3.d18*2000/22000**2) 
c            ppbest(4) = ppbest(4)+dlog10(dm)-dlog10(Lsol)           ! log Lk
            ppbest(2) = libppara(21,reclib)+dlog10(dm)              ! log Lnuv/Lo
            ppbest(3) = libppara(22,reclib)+dlog10(dm)              ! log Lr/Lo
            ppbest(4) = libppara(23,reclib)+dlog10(dm)              ! log Lk/Lo
c
            ppbest(5) = libppara(5,reclib)+dlog10(dm)               ! log Ldust/Lo
c
            ppbest(6) = 0.d0
            if (libppara(6,reclib).gt.0) 
     >       ppbest(6) = dlog10(libppara(6,reclib))+dlog10(dm)      ! Log Mass 

            ppbest(7) = -30.d0
            if (libppara(7,reclib).gt.0) 
     >      ppbest(7) = dlog10(libppara(7,reclib))+dlog10(dm)       ! Log SFR1e8 

            ppbest(8) = -30.d0
            if (libppara(8,reclib).gt.0) 
     >      ppbest(8) = dlog10(libppara(8,reclib))+dlog10(dm)       ! Log SFR1e9

            ppbest(9) = 0.d0
            if (libppara(9,reclib).gt.0) 
     >       ppbest(9) = dlog10(libppara(9,reclib))                 ! log gamma(Gyr)

            ppbest(10)= (libppara(10,reclib))                       ! metallicity
            ppbest(11)= (libppara(11,reclib))                       ! tauv 
            ppbest(12)= (libppara(12,reclib))                       ! mu

            ppbest(13) = 0.d0
            if (libppara(13,reclib).gt.0) 
     >      ppbest(13)= dlog10(libppara(13,reclib))                 ! log Age  [yr]  

            ppbest(14) = 0.d0
            if (libppara(14,reclib).gt.0) 
     >      ppbest(14)= dlog10(libppara(14,reclib))                 ! log Age  [yr] 

            ppbest(15) = 0.d0
            if (libppara(15,reclib).gt.0) 
     >      ppbest(15)= dlog10(libppara(15,reclib))                 ! log tlastburst [yr]  

            ppbest(16)= (libppara(16,reclib))                       ! frac burst 1.e8yr   
            ppbest(17)= (libppara(17,reclib))                       ! frac burst 1.e9yr 
            ppbest(18)= (libppara(18,reclib))                       ! Nb burst 

            ppbest(19) = 0.d0
            if (libppara(6,reclib).gt.0.and.libppara(7,reclib).ge.0)
     >      ppbest(19)= dlog10(libppara(7,reclib))-
     >                  dlog10(libppara(6,reclib))                  ! SFR(1e8)/Mass

            ppbest(20) = 0.d0
            if (libppara(6,reclib).gt.0.and.libppara(8,reclib).ge.0)
     >      ppbest(20)= dlog10(libppara(8,reclib))-
     >                  dlog10(libppara(6,reclib))                  ! SFR(1e9)/Mass

            ppbest(21) = libppara(24,reclib)                        ! A(FUV)  Mag 
            ppbest(22) = libppara(25,reclib)                        ! A(NUV)  Mag 
            ppbest(23) = libppara(26,reclib)                        ! A(G)    Mag 
            ppbest(24) = libppara(27,reclib)                        ! A(R)    Mag 
            ppbest(25) = libppara(28,reclib)                        ! A(Z)    Mag 
            ppbest(26) = libppara(29,reclib)                        ! A(K)    Mag 
            ppbest(27) = libppara(30,reclib)+dlog10(dm)             ! log Ldust/Lo
            ppbest(28) = libppara(31,reclib)                        ! D4000
         endif    

cccc    evolution through the dynamical range 
         if (chi2_phys.lt.1e5 .and. chi2_phys.gt. 0
     >      .and. libppara(1,reclib).gt.0 .and.npsum.gt.0) then
           
           do j = 1, 28
             do nchi = 1,lirmax(j)                    
c log  age
              if ( j.eq.1 .and. nppara(1).eq.1 .and.
     >    SNGL(DABS(dlog10(libppara(1,reclib))
     >            -physpmin(1,nchi))).le.SNGL(dbin(1)/2.d0) )
     >        chiphys(1,nchi)=chiphys(1,nchi)+dexp(-0.5*chi2_phys)
c log lum Lnuv
             if ( j.eq.2 .and. nppara(2).eq.1 .and.
     >    SNGL(DABS(libppara(21,reclib)+dlog10(dm)
c     >             -dlog10(Lsol*3.d18*400/2300**2)
     >             -physpmin(2,nchi))).le. SNGL(dbin(2)/2.d0) )
     >       chiphys(2,nchi)=chiphys(2,nchi)+dexp(-0.5*chi2_phys)
c log lum Lr
             if ( j.eq.3 .and. nppara(3).eq.1 .and.
     >    SNGL(DABS(libppara(22,reclib)+dlog10(dm)
c     >             -dlog10(Lsol*3.d18*1000/6000**2)
     >             -physpmin(3,nchi))).le.SNGL(dbin(3)/2.d0) )
     >       chiphys(3,nchi)=chiphys(3,nchi)+dexp(-0.5*chi2_phys)
c log lum Lk
             if ( j.eq.4 .and. nppara(4).eq.1 .and. 
     >    SNGL(DABS(libppara(23,reclib)+dlog10(dm)
c     >             -dlog10(Lsol*3.d18*2000/22000**2)
     >             -physpmin(4,nchi))).le.SNGL(dbin(4)/2.d0) )
     >       chiphys(4,nchi)=chiphys(4,nchi)+dexp(-0.5*chi2_phys)

c log ldust
             if ( j.eq.5 .and. nppara(5).eq.1 .and. 
     >    SNGL(DABS(libppara(5,reclib)+dlog10(dm)
     >             -physpmin(5,nchi))).le.SNGL(dbin(5)/2.d0) )
     >       chiphys(5,nchi)=chiphys(5,nchi)+dexp(-0.5*chi2_phys)

c log mass
             if ( j.eq.6 .and. nppara(6).eq.1 .and.
     >    SNGL(DABS(dlog10(libppara(6,reclib))+dlog10(dm)
     >             -physpmin(6,nchi))).le.SNGL(dbin(6)/2.d0) ) 
     >        chiphys(6,nchi)=chiphys(6,nchi)+dexp(-0.5*chi2_phys)
c log  sfr
             if ( j.eq.7 .and. nppara(7).eq.1
     >                   .and. libppara(7,reclib).gt.0 .AND.
     >    SNGL(DABS(dlog10(libppara(7,reclib)) +dlog10(dm)
     >             -physpmin(7,nchi))).le.SNGL(dbin(7)/2.d0) )
     >        chiphys(7,nchi)=chiphys(7,nchi)+dexp(-0.5*chi2_phys)
c  log sfr 
             if ( j.eq.8 .and. nppara(8).eq.1
     >                   .and. libppara(8,reclib).gt.0 .AND. 
     >    SNGL(DABS(dlog10(libppara(8,reclib))+dlog10(dm)
     >             -physpmin(8,nchi))).le.SNGL(dbin(8)/2.d0) )
     >        chiphys(8,nchi)=chiphys(8,nchi)+dexp(-0.5*chi2_phys)
c  log age
              if ( j.eq.9 .and. nppara(9).eq.1 .and. 
     >    SNGL(DABS(dlog10(libppara(9,reclib))
     >            -physpmin(9,nchi))).le.SNGL(dbin(9)/2.d0) )
     >        chiphys(9,nchi)=chiphys(9,nchi)+dexp(-0.5*chi2_phys)
c met
             if ( j.eq.10 .and. nppara(10).eq.1 .and.
     >    SNGL(DABS(libppara(10,reclib)
     >            -physpmin(10,nchi)))  .le. SNGL(dbin(10)/2.d0) )
     >       chiphys(10,nchi)=chiphys(10,nchi)+dexp(-0.5*chi2_phys)
c tauv
             if ( j.eq.11 .and. nppara(11).eq.1 .and.
     >    SNGL(DABS(libppara(11,reclib)
     >           -physpmin(11,nchi))) .le.SNGL(dbin(11)/2.d0) )
     >       chiphys(11,nchi)=chiphys(11,nchi)+dexp(-0.5*chi2_phys)
c mu
             if ( j.eq.12 .and. nppara(12).eq.1 .and.
     >    SNGL(DABS(libppara(12,reclib)
     >           -physpmin(12,nchi))).le.SNGL(dbin(12)/2.d0) )
     >       chiphys(12,nchi)=chiphys(12,nchi)+dexp(-0.5*chi2_phys)
c log  age
             if ( j.eq.13 .and. nppara(13).eq.1 .and.
     >    SNGL(DABS(dlog10(libppara(13,reclib))
     >            -physpmin(13,nchi))).le.SNGL(dbin(13)/2.d0) )
     >        chiphys(13,nchi)=chiphys(13,nchi)+dexp(-0.5*chi2_phys)
c log  age
              if ( j.eq.14 .and. nppara(14).eq.1 .and.
     >    SNGL(DABS(dlog10(libppara(14,reclib))
     >            -physpmin(14,nchi))).le.SNGL(dbin(14)/2.d0) )
     >        chiphys(14,nchi)=chiphys(14,nchi)+dexp(-0.5*chi2_phys)
c  log age
              if ( j.eq.15 .and. nppara(15).eq.1 .and.
     >    SNGL(DABS(dlog10(libppara(15,reclib))
     >            -physpmin(15,nchi))).le.SNGL(dbin(15)/2.d0) )
     >        chiphys(15,nchi)=chiphys(15,nchi)+dexp(-0.5*chi2_phys)
c  frac 1e8 
              if ( j.eq.16 .and. nppara(16).eq.1 .and.
     >    SNGL(DABS(libppara(16,reclib)
     >            -physpmin(16,nchi))).le.SNGL(dbin(16)/2.d0) )
     >        chiphys(16,nchi)=chiphys(16,nchi)+dexp(-0.5*chi2_phys)
c  frac 1e9
              if ( j.eq.17 .and. nppara(17).eq.1 .and.
     >    SNGL(DABS(libppara(17,reclib)
     >            -physpmin(17,nchi))).le.SNGL(dbin(17)/2.d0) )
     >        chiphys(17,nchi)=chiphys(17,nchi)+dexp(-0.5*chi2_phys)
c  Nb burst 
              if ( j.eq.18 .and. nppara(18).eq.1 .and.
     >    SNGL(DABS(libppara(18,reclib)
     >            -physpmin(18,nchi))) .le.SNGL(dbin(18)/2.d0) )
     >        chiphys(18,nchi)=chiphys(18,nchi)+dexp(-0.5*chi2_phys)
c log SSFR1e8
             if ( j.eq.19 .and. nppara(19).eq.1 .and.
     >    SNGL(DABS(dlog10(libppara(7,reclib))
     >            -dlog10(libppara(6,reclib))
     >            -physpmin(19,nchi))).le.SNGL(dbin(19)/2.d0) )
     >        chiphys(19,nchi)=chiphys(19,nchi)+dexp(-0.5*chi2_phys)
c  log SSFR1e9
             if ( j.eq.20 .and. nppara(20).eq.1 .and. 
     >    SNGL(DABS(dlog10(libppara(8,reclib))
     >           -dlog10(libppara(6,reclib))
     >           -physpmin(20,nchi))).le.SNGL(dbin(20)/2.d0) )
     >        chiphys(20,nchi)=chiphys(20,nchi)+dexp(-0.5*chi2_phys)
c  A fuv 
             if ( j.eq.21 .and. nppara(21).eq.1 .and.
     >    SNGL(DABS(libppara(24,reclib)
     >           -physpmin(21,nchi))) .le.SNGL(dbin(21)/2.d0) )
     >       chiphys(21,nchi)=chiphys(21,nchi)+dexp(-0.5*chi2_phys)
c  A nuv 
             if ( j.eq.22 .and. nppara(22).eq.1 .and.
     >    SNGL(DABS(libppara(25,reclib)
     >           -physpmin(22,nchi))) .le.SNGL(dbin(22)/2.d0) )
     >       chiphys(22,nchi)=chiphys(22,nchi)+dexp(-0.5*chi2_phys)
c  A g 
             if ( j.eq.23 .and. nppara(23).eq.1 .and.
     >    SNGL(DABS(libppara(26,reclib)
     >           -physpmin(23,nchi))) .le.SNGL(dbin(23)/2.d0) )
     >       chiphys(23,nchi)=chiphys(23,nchi)+dexp(-0.5*chi2_phys)
c  A r 
             if ( j.eq.24 .and. nppara(24).eq.1 .and.
     >    SNGL(DABS(libppara(27,reclib)
     >           -physpmin(24,nchi))) .le.SNGL(dbin(24)/2.d0) )
     >       chiphys(24,nchi)=chiphys(24,nchi)+dexp(-0.5*chi2_phys)
c  A z 
             if ( j.eq.25 .and. nppara(25).eq.1 .and.
     >    SNGL(DABS(libppara(28,reclib)
     >           -physpmin(25,nchi))) .le.SNGL(dbin(25)/2.d0) )
     >       chiphys(25,nchi)=chiphys(25,nchi)+dexp(-0.5*chi2_phys)
c  A k 
             if ( j.eq.26 .and. nppara(26).eq.1 .and.
     >    SNGL(DABS(libppara(29,reclib)
     >           -physpmin(26,nchi))) .le.SNGL(dbin(26)/2.d0) )
     >       chiphys(26,nchi)=chiphys(26,nchi)+dexp(-0.5*chi2_phys)
c  Ldust Hot 
             if ( j.eq.27 .and. nppara(27).eq.1 .and. 
     >    SNGL(DABS(libppara(30,reclib)+dlog10(dm)
     >             -physpmin(27,nchi))).le.SNGL(dbin(27)/2.d0) )
     >       chiphys(27,nchi)=chiphys(27,nchi)+dexp(-0.5*chi2_phys)
c  D4000 
             if ( j.eq.28 .and. nppara(28).eq.1 .and. 
     >    SNGL(DABS(libppara(31,reclib)
     >             -physpmin(28,nchi))).le.SNGL(dbin(28)/2.d0) )
     >       chiphys(28,nchi)=chiphys(28,nchi)+dexp(-0.5*chi2_phys)
c
           enddo
          enddo
        endif 
 6      continue
      enddo
 4    close(29)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   measuring the median and quartiles 
c      write(UO,*) ' analysis 0'
      valtest=-99    
      do k = 1,28
        ppmed(k) = -99.
        ppinf(k) = -99.   
        ppsup(k) = -99.  
        if (nppara(k).eq.1 .and. chirmin.le.10000.d0 ) then 
          test=0
          do j = 1,lirmax(k)
            if (chiphys(k,j).lt.1.d-60) chiphys(k,j)=0.d0
            lumir(j)=physpmin(k,j)
            chilir(j)=chiphys(k,j)
            if (chilir(j).gt.0) test = 1
          enddo
          if (test.gt.0) then 
            bstart= sbin(k)
            bend  = ebin(k)
            lirmaxs=lirmax(k)
            call TRAPZD(lumir,chilir,lirmaxs,bstart,bend,barea)
            if (barea.gt.0) then 
              call PROB_MED(lumir,chilir,lirmaxs,barea,zmed,zpdzi,zpdzs)
               ppmed(k) = zmed
               ppinf(k) =zpdzi   
               ppsup(k) =zpdzs     
            endif
          endif
          valtest=ppmed(6)
        endif
      enddo 
c      write(*,*) ' '
c      write(*,*) chirmin,reclphysb
c      write(*,*) ' ' 
cccccccccccccccccccccccccccccccc
c  test on one of the parameter to see if measured 
c  CANNOT AFFORD FOR A SECOND LOOP ON 100000 models ! 
c      if ( iter.eq. 0 .AND. nf_phys.gt.3 .AND. 
c     >   ((chirmin/(nf_phys-3)).gt.5 .OR. valtest.eq.-99.) ) then
      if (iter.eq.0 .and. nf_phys.gt.3 .AND.
     >      (valtest.lt.-99 .OR. chirmin.gt.250) ) then 
         iter=iter+1
         facerr=dsqrt(chirmin/(nf_phys-3))
         write(UO,*) 'redoing with chi2=',chirmin,nf_phys,facerr
         do k = 1,imag
            if (iter.eq.1) sabo(k) = sab(k)   
            sab(k) = sab(k)*facerr
         enddo
         chiorig=chirmin 
         goto 3
      endif
c
      if ( iter.gt.0) then 
        do k = 1, imag
           sab(k) = sabo(k)
        enddo 
        chipbest=chiorig
c         write(UO,*) 'redone with chi2=',chirmin
      endif    
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c   extracting  mag at z=0 for the best model
      if (reclphysb.gt.0) then  
        open(29,file=file(1:lnblnk(file)),status='unknown',
     >       access='direct',recl=reclmax)
        do irec = zrecpi(1),zrecps(1)
           read(29,rec=irec) modlib,vali,val,zlib,
     >        val2,reclib,imag,(maglib(k),k=1,imag),
     >       (kcorlib(k),k=1,imag)
           if (modlib.eq.reclphysb-1) then 
             do k = 1, imag 
               magphys0(k) = maglib(k) -2.5*dlog10(dmpbest)
             enddo 
c             write(UO,'(f6.3,1x,2(i8,1x),2(f8.3,1x))')
c     >         zlib,modlib,reclphysb,magphys0(1),magphys0(3) 
             goto 5 
           endif
         enddo
 5       close(29)
      endif     
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
 2    return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Return zmin,zmax included inside 50%,68%,90%,95% of the area 
c
c      PROGRAM PROBAZBAY
      SUBROUTINE PROB_MED(x,y,bmax,area,zmed,zpdzi,zpdzs)
      implicit none 
c      
      INTEGER*4 bmax,i
      REAL*8    x(500),y(500),yn(500),zmed,slope
      REAL*8    zli,zls,area,proba,probai(500),zpdzi,zpdzs
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
