c       
c       program    SC_LIB
c
	subroutine READ_SCLIB(imodmax,magmod,paramod,parasc_lim,izmax) 
c              
c   * Reads the S Charlot MODELs at various Z for fixed bands
c     as defined in SC models (BVRIJK)VVDS + (ugriz)CFHTLS + (IRAC1234)    
c     build the list of models Mag according to filters :
c     FILTER_FILE : /data/arnouts/zpwork/filt/vvds_cfhtls_swire.filt
c     U B V R I up gp rp ip zp J Ks IRAC_1 IRAC_2 IRAC_3 IRAC_4 24mic
c      Output format 
c      Model Z  Age Age_wr TauV mu  M*  M*/(1-R) SFR(1e8)
c      Return fluxes without the Distance modulus
c      U B V R I u g r i z J K IRAC1 2 3 4 24mic 
        implicit none 
c   Array declarations
      character*512   name,scmod
c
      integer*4   test,lnblnk
      integer*4   k,isc,iz,izmax
      integer*4   nmod
      character*512  zpwork,zpdir
      character*1024 str
c
      real*8    para(10),magm(15)
      integer*4   sclib,iscmax
c  all in one file
c      PARAMETER   (sclib=800000)
c      real*8    paramod(10,sclib),parasc_lim(2,10)     
c      real*8    zmod(sclib),magmod(17,sclib) 
c  one file per z 
      PARAMETER   (sclib=90000)
      real*8    paramod(15,10,sclib),parasc_lim(2,10)     
c      real*8    zmod(15,sclib),magmod(15,17,sclib) 
      real*8    magmod(15,17,sclib) 
      integer*4 imodmax(15)
c
      real*8    funz,dm,z,h0,om0,l0
      external funz

      write(*,*) ''
      write(*,*) ''
c   cosmology   
      h0=70
      om0=0.3
      l0=0.7 
c
cccccccccccccccccccccccccccccccccccccccccccc
c   environmental variable 
      call getenv('LEPHAREDIR',zpdir)
      test=lnblnk(zpdir)
      if (test .eq. 0) then
        write(6,*) 'ERROR : variable LEPHAREDIR not defined --> STOP'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(6,*) 'ERROR : variable LEPHAREWORK not defined --> STOP'
        stop
      endif
cccccccccccccccccccccccccccccccccccccccccccc
c   reading SC models 
      do k = 1, 10
         parasc_lim(1,k) = 1.e10
         parasc_lim(2,k) = -1.e10
      enddo
      scmod = zpdir(1:lnblnk(zpdir))//'/sed/GAL/SC_MODELS/SC_MOD.list'
      open(1,file=scmod,status='old',err=41)
      nmod = 0
      isc=0
      iz = 0
      izmax = 0
      do while (.true.)
         read(1,'(A)',end=2) str
         if (str(1:1).ne.'#' .AND. str(1:1).ne.' ' 
     >      .AND. str(1:1).ne.char(9)
     >      .AND. str(1:1).ne.char(13)) then
	    nmod = nmod + 1 
            z  = 0.1*float(nmod)
            dm = funz(z,h0,om0,l0)
c
c
            name=zpdir(1:lnblnk(zpdir)) // '/sed/GAL/SC_MODELS/'
     >                 // str(1:lnblnk(str))
cccccccccccccccc
c   if ASCII file (A)
	    open (21,file=name,status='old',err=42)
            read(21,*)
            read(21,*)  
            iz = 0 
            do while (.true.)
c            read(21,*,end=31) age,ager,tauv,mu,mstar,mstarr,sfr,(magm(k),k=1,15) 
              read(21,*,end=31) (para(k),k=1,7),(magm(k),k=1,15) 
              isc = isc + 1
              iz = iz + 1
c              zmod(isc)=z
c              zmod(nmod,iz) = z
              do k = 1, 7
                if (k.ge.3 .and. k.le.4) then 
c                   paramod(k,isc) = para(k)
                   paramod(nmod,k,iz) = para(k)                   
                else
c                   paramod(k,isc) = dlog10(para(k))
                   paramod(nmod,k,iz) = dlog10(para(k))
                endif
c             if (paramod(k,isc).le.parasc_lim(1,k))
c     >                            parasc_lim(1,k)=paramod(k,isc)
c             if (paramod(k,isc).ge.parasc_lim(2,k))
c     >                            parasc_lim(2,k)=paramod(k,isc)
             if (paramod(nmod,k,iz).le.parasc_lim(1,k))
     >                            parasc_lim(1,k)=paramod(nmod,k,iz)
             if (paramod(nmod,k,iz).ge.parasc_lim(2,k))
     >                            parasc_lim(2,k)=paramod(nmod,k,iz)
              enddo
c             (UBVRI) VVDS     
              magmod(nmod,1,iz) = -99
              magmod(nmod,2,iz) = magm(1)  ! + dm
              magmod(nmod,3,iz) = magm(2)  ! + dm
              magmod(nmod,4,iz) = magm(3)  ! + dm 
              magmod(nmod,5,iz) = magm(4)  ! + dm 
c            (ugriz) CFHTLS
              magmod(nmod,6,iz) = magm(7)  ! + dm
              magmod(nmod,7,iz) = magm(8)  ! + dm
              magmod(nmod,8,iz) = magm(9)  ! + dm 
              magmod(nmod,9,iz) = magm(10) ! + dm 
              magmod(nmod,10,iz)= magm(11) ! + dm 
c            (JK) VVDS
              magmod(nmod,11,iz) = magm(5) ! + dm
              magmod(nmod,12,iz) = magm(6) ! + dm

c            (IRAC 1 2 3 4 ) SWIRE 
              magmod(nmod,13,iz) = magm(12) ! + dm
              magmod(nmod,14,iz) = magm(13) ! + dm
              magmod(nmod,15,iz) = magm(14) ! + dm
              magmod(nmod,16,iz) = magm(15) ! + dm

c            (24mic) SWIRE
              magmod(nmod,17,iz) = -99
c     converting mag to flux in AB 
              do k = 2,16 
                if (magmod(nmod,k,iz) .gt. 99) then
                 magmod(nmod,k,iz) = 0.
                else
                 magmod(nmod,k,iz)=10**(-0.4*(magmod(nmod,k,iz)+48.59)) 
                endif 
              enddo
c              if (mod(isc,30000).eq.0) then 
c               write(*,'(100(E12.6,1x))') (magmod(k,isc),k=1,17)
c             write(*,'(10(f8.3,1x))') zmod(isc),(paramod(k,isc),k=1,7)
c              endif
            enddo
 31	    close(21)
            imodmax(nmod) = iz 
            if (iz.gt.izmax) izmax = iz 
         write(*,'(A,2x,f3.1,2x,f5.2,1x,i8)') str(1:lnblnk(str)),z,dm,iz
	 endif  
      enddo
 2    close(1) 
      iscmax = isc  
c    print min max for parameters 
       do k = 1,7
          write(*,*) parasc_lim(1,k),parasc_lim(2,k)
       enddo
c
        write(*,'(A,1x,i8)') 'izmax = ',izmax
c
        goto 20
 41	write (6,*) 'File ',scmod(1:lnblnk(scmod)),' not found'
 42	write (6,*) 'File ',name(1:lnblnk(name)),' not found'
c
 20	return 
        end
c
