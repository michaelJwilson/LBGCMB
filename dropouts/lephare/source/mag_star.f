c   last modif : 19/10/00
      program mag_star
c
c     Construction des differents modeles d'etoiles pour la 
c     mesure des couleurs d'etoiles et pour etre utiliser dans 
c     la mesure des z-photometriques 
c 
c      Librairie stellaire vient de Pickles, 1998, PASP 110 
c
      implicit none
      integer*4 nbf,maxsize
      integer*4 wmax
      parameter (maxsize=110000)
c      parameter (wmax=8000)
      INCLUDE 'dim_wave.decl'      
      INCLUDE 'dim_filt.decl'
      INCLUDE 'out_unit.decl' 
      integer*4 jmax(nbf),id(nbf)
      integer*4 imag,imas,i,j,k,if,iw,nlrep
      integer*4 recmax,reclmax,irec
      integer*4 nmod,nr,orec,orecmax
      integer*4 test,dummy1,pass
      integer*4 nlveg,lnblnk
c
      real*8    wave(wmax),dummy,vec(wmax)
      real*8    dummy2
      real*8    magveg,magab,abveg,mag(nbf),abcor(nbf)
      real*8    flux(wmax)
      real*8    lambf(nbf,maxsize),repf(nbf,maxsize)
      real*8    lamb(maxsize),rep(maxsize)
      real*8    lbveg(wmax),flveg(wmax)
      real*8    zp(nbf),flmoy(nbf),flwidth(nbf)
      real*8    fcorr(nbf)

c
      character*4096  magtyp,name(nbf),param
      character*4096  paravc(500)
      character*4096  vega,valc
      character*4096  str
      character*4096  file,filters,lib,colib,config,outasc
      character*4096  zpdir,zpwork
      character*4096  cr_date,fdate

c      If you want to give a name to the output screen file
      if(UO.eq.30)
     .      open(UO,file='/tmp/screenMagStar.dat',status='unknown') 
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c  environmental variable 
      call getenv('LEPHAREDIR',zpdir)
      test=lnblnk(zpdir)
      if (test .eq. 0) then
        write(UO,*) 'WARNING :  variable LEPHAREDIR not defined'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(UO,*) 'WARNING :  variable LEPHAREWORK not defined'
        stop
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c help on line
      param='mag_star'
      call get_help(test)
      if (test .eq. 1) call help(param)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c Initialisation of Input parameter from config file
c  read option
      param='-c'
      call get_conf(param,config,test)
      if (test.ne.1) call err_option(param,1)
      call get_path(config)
      param='-FILTER_FILE'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  filters=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
      param='-MAGTYPE'
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1) magtyp=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
      param='-STAR_LIB_IN'
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1)  lib=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
      param='-STAR_LIB_OUT'
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1)  colib=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
      param='-LIB_ASCII'
      call  getc_option(param,config,1,paravc,test)
      if (test.eq.1)  outasc=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,2)
      if (test.ne.1)  outasc='NO'
cccccccccccccccccccccccccccc
c  read Vega spectrum
      vega = zpdir(1:lnblnk(zpdir)) // '/vega/VegaLCB.sed'
c       vega = zpdir(1:lnblnk(zpdir)) // '/vega/a0v_n.sed'
        open(1, file=vega,status='old')
        read(1,*)   
        i=0
	do while (.true.)
	  i=i+1
	  read(1,*,end=10)  lbveg(i), flveg(i)
 	enddo
 10     nlveg=i-1  
        close(1)
c  read the filters from  file
      file = zpwork(1:lnblnk(zpwork)) //'/filt/'
     >        // filters(1:lnblnk(filters))
      open(1,file=file,status='unknown',err=56)
      read(1,'(A)') str 
      call val_string(str,paravc,test)
      read(paravc(2),'(i10)') imag
      write(*,*) imag
      do i =  1, imag
         read(1,'(A)') str
         call val_string(str,paravc,test)
         read(paravc(2),'(i10)') jmax(i)
         read(paravc(5),'(i10)') id(i)
         name(i)=paravc(3)
         do j = 1,jmax(i)
	    read(1,*)  lambf(i,j), repf(i,j)
 	 enddo
      enddo
      close(1) 
      call zeropoint(filters,zp,abcor,flmoy,flwidth,fcorr,imag)
c
cccccccccccccccccccccccccccc
c  writing options 
       write(UO,'(A)')"###########################################"
       write(UO,'(A)')"# It s computing the SYNTHETIC MAGNITUDES #"
       write(UO,'(A)')"# For  STARS with the following OPTIONS   #"
       file = zpwork(1:lnblnk(zpwork)) //'/filt/'
     >      // filters(1:lnblnk(filters))
       write(UO,'(2A)')"# FILTER_FILE  : ",file(1:lnblnk(file))  
       write(UO,'(2A)')"# MAGTYPE      : ",magtyp(1:lnblnk(magtyp))
       write(UO,'(A,50(f8.3,1x))')
     >                 "# FLUX_COR    : ",(fcorr(j),j=1,imag)
        file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >      // lib(1:lnblnk(lib)) //'(.doc & .bin)'
       write(UO,'(2A)')"# STAR_LIB_IN  : ",file(1:lnblnk(file))
        file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >      // colib(1:lnblnk(colib)) //'(.doc & .bin)'
       write(UO,'(2A)')"# STAR_LIB_OUT : ",file(1:lnblnk(file))
       write(UO,'(2A)')"# LIB_ASCII    : ",outasc(1:lnblnk(outasc))  
       write(UO,'(A)')"###########################################"
cccccccccccccccccccccccccccc
c      
c  read number of CSP model and number of records per model nrec(i)
      file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >       // lib(1:lnblnk(lib)) //'.doc'
      param='NUMBER_SED'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') imas
      param='NUMBER_ROWS'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') recmax
      param='RECORD_LENGTH'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') reclmax
c  
      write(UO,*) 'number of record :',imas,recmax,reclmax
c
      file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >       // lib(1:lnblnk(lib)) //'.bin'
c  read wavelengthes and fluxes for each model 
c  irec    = no records for input file (1) 'lib'
c  orec    = no records for output file (2) 'colib'
c  reclmax = record length in input library
      open(1,file=file,status='unknown',access='direct',recl=reclmax)
c  Open output binary file with magnitude
      orecmax = 8*(3+imag)
      file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >      // colib(1:lnblnk(colib)) // '.bin'
      open(2,file=file,form='unformatted',status='unknown',
     > access='direct',recl=orecmax)
c  Open doc file with info for the binary file
      file = zpwork(1:lnblnk(zpwork)) //'/lib_mag/'
     >       // colib(1:lnblnk(colib)) // '.doc'
      open(4,file=file,status='unknown')
c  Open ASCII file with magnitude
      if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y") then
        file = colib(1:lnblnk(colib)) // '.dat'
        open(5,file=file,status='unknown')
        write(5,1005) '#',imag,magtyp(1:4),(name(j),j=1,imag)
      endif
      write(UO,1005) '*',imag,magtyp(1:4),(name(j),j=1,imag)
c  used for consistance with galaxy and QSO libraries
      dummy1= 1
      dummy2= 0.
      irec = 0
      orec = 0
      pass = 0 
c  computing magnitudes for each SED spectrum in all the filters
      do i = 1, recmax 
         irec = irec + 1
         read(1,rec=irec) nr,nmod,dummy,iw,(vec(k),k=1,iw)
         if (dummy.eq.-1.0) then
c  read lambda SEDs
c            write(6,*) irec,nr,nmod,dummy,iw,' lambda vector scaling'
            do k = 1,iw
               wave(k) = vec(k)
            enddo
         else
            pass = pass + 1
            write(UO,1010) nmod,irec,nr,dummy,iw,char(13)
            call flush(6)
c  read flux SEDs
            do k = 1, iw
               if (vec(k).lt.0) then
                  write(UO,1011) nmod,i,wave(k),vec(k),char(13)
                  call flush(6)
                  flux(k) = 0
               else
c    Flux defined as fl = F/Fvega(5556) in Pickles 1998
c    with Fvega(5560)=3.432E-09 erg.cm-2.s-1.A-1 done in sedtolib.f
                  flux(k) = vec(k)
               endif 
            enddo
c  compute mag for each filter 
           do j = 1, imag 
              nlrep = jmax(j)
              do k = 1,nlrep
	        lamb(k) = lambf(j,k)
                rep(k) =  repf(j,k)
 	      enddo
c   write : AB correction during the first pass
             if (pass.eq.1) then
               valc='V'
               call cal_mag(flux,wave,iw,rep,lamb,nlrep,
     >         flveg,lbveg,nlveg,valc,magveg,magab,abveg)
               abcor(j) = abveg            
             endif
c   computes the magnitude
             call cal_mag(flux,wave,iw,rep,lamb,nlrep,
     >          flveg,lbveg,nlveg,magtyp,magveg,magab,abveg)
c    check if lambda star pass fully through the filter
             if (lamb(nlrep).le.wave(iw)) then
                if (magtyp(1:1) .eq. 'A')   mag(j)= magab 
                if (magtyp(1:1) .eq. 'V')   mag(j)= magveg 
              else
                 mag(j)=999
              endif
c     apply CALIB correction for long wavelengths
              mag(j) = mag(j) -2.5*dlog10(fcorr(j))
           enddo
           if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y") then
              if (pass.eq.1) then
                write(5,1006) '#',imag,magtyp(1:4),
     >                 (abcor(if),if=1,imag) 
              endif     
c   write :  nmod,iext,zoss,age,imag,(mag(if),if=1,imag)
             write(5,1001) nmod,nr,imag,(mag(if),if=1,imag)
           endif
           orec = orec + 1
           write(2,rec=orec) nmod,nr,imag,(mag(if),if=1,imag)
         endif
      enddo
      write(UO,*) "                     "
      write(UO,*) " DONE       "
c
      cr_date=fdate()
      write(4,'(2A)')     "CONFIG_FILE      ",config(1:lnblnk(config))
      write(4,'(2A)')     "LIB_TYPE         ","STAR"
      write(4,'(2A)')     "LIB_NAME         ",lib(1:lnblnk(lib))
      write(4,'(A,3x,I8)')"NUMBER_SED       ",imas
      write(4,'(A,3x,I8)')"NUMBER_ROWS      ",orec
      write(4,'(A,3x,I8)')"RECORD_LENGTH    ",orecmax      
      write(4,'(51(A,2x))')"FILTER_FILE      ",
     >    filters(1:lnblnk(filters))
      write(4,'(51(A,2x))')    "FILTERS        ",
     >                           (name(j)(1:lnblnk(name(j))),j=1,imag)
      write(4,'(51A)')    "MAG_TYPE         ",magtyp(1:5)
      write(4,'(A,100(f8.3,1x))')
     >                    "AB_COR           ",(abcor(j),j=1,imag)
      write(4,'(A,50(f8.3,1x))')
     >                     "FLUX_COR         ",(fcorr(j),j=1,imag)
      write(4,'(2A)')     "CREATION_DATE    ",cr_date(1:lnblnk(cr_date))
      close(1)
      close(2)
      close(4)
      if (outasc(1:1).eq."Y" .or. outasc(1:1).eq."y") close(5)
      if(UO.eq.30) close(UO)
 1001 format(3(i6,1x),1000(f9.4,1x))
c 1002 format(I3,1x,A5,1x,60(A10,1x))
c 1003 format(I3,1x,A5,1x,60(f6.3,1x))
c 1004 format(A25,1x,A20)
 1005 format(A1,1x,I3,1x,A5,1x,100(A10,1x))
 1006 format(A1,1x,I3,1x,A5,1x,60(f6.3,1x))
 1010 format("model -> ",I6,1x,I8,1x,I6,1x,f8.3,1x,I6,1x,a1,$)
 1011 format("model -> ",I6,1x,I8,1x,"Neg. flux in ",
     >       f8.4,1x,E12.6,1x,a1,$)
c
      stop
 56   write (6,*) 'File ',file(1:lnblnk(file)),
     > ' not found -> STOP '

      end
c
