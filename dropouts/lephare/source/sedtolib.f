c       
        PROGRAM sedtolib 
c   2 options : sedtolib 
c                        -c zphot.para          : config file 
c                        -t [G/g or Q/q or S/s] : Gal/QSO/STAR list
c              
c   * Reads the SED  list for Gal and QSO and STARS  
c     from ASCII files     : A for STARS/QSO/GAL with no-evolution
c                            or F (ASCII from PEGASE spectra)
c     or from BINARY files : B for GISSEL spectra 
c   * Writes the SEDs  in BINARY LIBRARY
c     and writes the associated DOC files (in ASCII format)
c    
c   We uses only the number of ages given in a file (SEL_AGE parameter) 
c
c   WARNING : Maximum number of wavelength elements : 8000 
c  
        implicit none 
 
        INCLUDE 'out_unit.decl' 

c   Array declarations
        character*4096   name,mod,lib,doc,name2
c   Maximum number of wavelength points 
        integer*4   imw,wmax,jtype
        integer*4   test,lnblnk
c	PARAMETER   (wmax=8000)
        INCLUDE 'dim_wave.decl'      
        integer*4   i,j,k,n,nsteps,inw,kmax
        integer     nmod,reclmax,recpmax,nrec,nrec0
        integer     cage(0:500),nage     ! ,nobj
	real*4      hgs(wmax),tbs(0:500),ws(wmax)
	real*8      h(wmax),fluxconv
	real*8      hg(wmax),hg0(wmax)
	real*8      w(wmax),tb(0:500),dummy
        real*8      newage(500),paravr(500),fscale
        real*8      ltirg(500),luvg(500),loptg(500)
        real*8      lnirg(500),d4000g(500)
        real*8      agemin,agemax
        real*8      hir(wmax),wir(wmax),hir0(wmax)
        real*8      age,mass,sfr,tau,d4000
        real*8      ltir,ltir0,ltir1,lnir,lnirk,lfuv,lnuv,luv,lopt,loptr
        real*8      ldustc,ldusth,afuv,anuv,aopt,anir,a_g,a_z
        character*4096  zpdir,zpwork,config,param
        character*4096  fileage,typ,sedtyp,path,phys,libphys
	character*4096  paravc(500),str2,str
        character*4096  cr_date,fdate
c  variable from BC03
        integer*4   ix
        real*4      fx(100),bol(0:500),mstr(0:500),sf(0:500),evf(0:500)
        real*4      snr(0:500),pnr(0:500),bh(0:500),sn(0:500),wd(0:500)
        real*4      rm(0:500)
	logical     stelib
        integer   iseg,iop,imf
        parameter (imf=10)
        real        xx(imf),lm(imf),um(imf),baux(imf),cc(imf),cn(imf)
        real        ml,mu0,totm,totn,avs,jo,tau0,tau1,tau2,tau3,tau4
	character   id*80,id2*80,id3*80
c  variable for BC07_STOCHASTIC 
        REAL      tform,gamma,zmet,tauv0,mu,mstr1,mstr0
        REAL      tlastburst,fburst(5),ftot(5),age_wm,age_wr
        REAL      ages(900),sfrs(900),sfrav(900)
        REAL      fprop(wmax),fprop0(wmax)
        INTEGER   nburst,nsage,jstoch            ! ,nsagemax
        real      eta 

c  variable from PEGASE2 (FIOC et al) 
        integer*4   inlines
	real*8      wline(wmax)
        real*8      hline(wmax)
        real*8      ageref(500)
        integer*4   amax,found

c
	real*8      Lsol,pc,c
	parameter  (Lsol=3.826e33,pc=3.086e18)
        parameter  (c=2.99792458e+18)      

c      If you want to give a name to the output screen file
      if(UO.eq.30)
     .      open(UO,file='/tmp/screenSedtolib.dat',status='unknown') 

   
c
cccccccccccccccccccccccccccccccccccccccccccc
c   environmental variable 
      call getenv('LEPHAREDIR',zpdir)
      test=lnblnk(zpdir)
      if (test .eq. 0) then
        write(UO,*) 'ERROR : variable LEPHAREDIR not defined --> STOP'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(UO,*) 'ERROR : variable LEPHAREWORK not defined --> STOP'
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c help on line
      param='sedtolib'
      call get_help(test)
      if (test .eq. 1) call help(param)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Initialisation of Input parameter from config file
      param='-c'
      call get_conf(param,config,test)
      if (test.ne.1) call err_option(param,1)
      call get_path(config)
cccccccc  OPTIONS  ccccccccccccccccc
      param='-t'
      call get_conf(param,sedtyp,test)
      if (test.ne.1) call err_option(param,1)
c
      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
        sedtyp = 'G'
	path='/sed/GAL/'
c
        param='-GAL_SED'
        call getc_option(param,config,1,paravc,test)
        if (test.eq.1) mod = paravc(1)(1:lnblnk(paravc(1)))
        if (test.ne.1) call err_option(param,1)
        call get_path(mod)
c
        param='-GAL_FSCALE'
        call getf_option(param,config,1,paravr,test)
        if (test.eq.1) fscale = paravr(1)
        if (test.ne.1) call err_option(param,2)
        if (test.ne.1) then
          fscale=1
	  write(UO,'(A)') "option -GAL_FSCALE sets to 1"
        endif
c
        param='-GAL_LIB'
        call getc_option(param,config,1,paravc,test)
        if (test.eq.1) lib = paravc(1)(1:lnblnk(paravc(1)))
        if (test.ne.1) call err_option(param,1)
c
        param='-SEL_AGE'
        call getc_option(param,config,1,paravc,test)
	if (test.eq.1) fileage = paravc(1)(1:lnblnk(paravc(1)))
        if (test.ne.1) fileage = 'NONE'
        call get_path(fileage)
c
        param='-AGE_RANGE'
        call getf_option(param,config,2,paravr,test)
        if (test.eq.2) agemin=paravr(1)
        if (test.eq.2) agemax=paravr(2) 
        if (test.ne.2) agemin=0.
        if (test.ne.2) agemax=1.e12
c
      elseif (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q") then
        sedtyp = 'Q'
	path='/sed/QSO/'
c
        param='-QSO_SED'
        call getc_option(param,config,1,paravc,test)
        if (test.eq.1) mod = paravc(1)(1:lnblnk(paravc(1)))
        if (test.ne.1) call err_option(param,1)
        call get_path(mod)
c
        param='-QSO_FSCALE'
        call getf_option(param,config,1,paravr,test)
        if (test.eq.1) fscale = paravr(1)
        if (test.ne.1) call err_option(param,2)
        if (test.ne.1) then 
 	  fscale=1
	  write(UO,*) "option -QSO_FSCALE sets to 1"
        endif
c
        param='-QSO_LIB'
        call getc_option(param,config,1,paravc,test)
        if (test.eq.1) lib = paravc(1)(1:lnblnk(paravc(1)))
        if (test.ne.1) call err_option(param,1)
	fileage= ' ' 
c
        libphys = 'NO' 
c
      elseif (sedtyp(1:1).eq."s" .or. sedtyp(1:1).eq."S") then
        sedtyp = 'S'
	path='/sed/STAR/'
c
        param='-STAR_SED'
        call getc_option(param,config,1,paravc,test)
        if (test.eq.1) mod = paravc(1)(1:lnblnk(paravc(1)))
        if (test.ne.1) call err_option(param,1)
        call get_path(mod)
c
	param='-STAR_LIB'
        call getc_option(param,config,1,paravc,test)
        if (test.eq.1) lib = paravc(1)(1:lnblnk(paravc(1)))
        if (test.ne.1) call err_option(param,1)
	fileage= ' '
c
        param='-STAR_FSCALE'
        call getf_option(param,config,1,paravr,test)
        if (test.eq.1) fscale = paravr(1)
        if (test.ne.1) call err_option(param,2)
        if (test.ne.1) then 
 	  fscale=1
	  write(UO,*) "option -STAR_FSCALE sets to 1"
        endif
c
        libphys = 'NO' 
      else 
	  write(UO,*) 'Wrong setting for -t option'
	  write(UO,*) ' must be : -t G(g)  or -t Q(q) or -t S(s) '
          STOP
      endif
c 
cccccccccccccccccccccccccccccccccccccccccc  
c   Read the number of seds in file  mod and check the number of lambdas 
      write(UO,*) '   '
      write(UO,*) ' first pass : reading each SED ...'
c
      mod =  mod(1:lnblnk(mod))
      open(1,file=mod,status='old',err=41)
      nmod = 0
      imw  = 0
      jstoch = 0             ! number of models in stochastic library 
      do while (.true.)
         read(1,'(A)',end=2) str
         if (str(1:1).ne.'#' .AND. str(1:1).ne.' ' 
     >        .AND. str(1:1).ne.char(9)
     >        .AND. str(1:1).ne.char(13)) then
c
	    call val_string(str,paravc,test)
	    name=paravc(1)
            if (test .eq. 1) then 
              typ = 'A'
            else
              typ=paravc(2)
            endif
	    nmod = nmod + 1 
            name2=zpdir(1:lnblnk(zpdir)) // path(1:lnblnk(path))
     >                 // name(1:lnblnk(name))
cccccccccccccccc
c   if ASCII file (A) or if library with Long Wavelength  (LW)
	    if (typ(1:1) .eq. 'A' .or. typ(1:2).eq.'LW' ) then
	       open (21,file=name2,status='old',err=4)
	       i = 0
               do while (.true.)
                     i = i + 1
 	            read(21,*,end=31)  w(i),h(i)                  
               enddo
 31	       inw = i - 1
               close(21)
	       if (inw.gt.imw) imw=inw
cccccccccccccccc
c   if GISSEL file (B)
            elseif (typ(1:lnblnk(typ)) .eq. "B" .or.
     &               typ(1:lnblnk(typ)) .eq. "BC03") then
c               write(6,*) ' opening ',name2(1:lnblnk(name2))
	       open (21,file=name2,form='unformatted',status='old',err=4)
	       read (21)  nsteps,(tbs(i),i=0,nsteps-1)       
	       read (21)  inw,(ws(i),i=1,inw)
               close(21)
	       if (inw.gt.imw) imw=inw
cccccccccccccccc
c   if PEGASE file (F)
	    elseif (typ(1:1) .eq. "F") then 
 	       open (21,file=name2,status='old',err=4)
               read(21,'(a)') str  
               do while (str(1:10).ne.'**********')
                 read(21,'(a)') str 
               end do
               read(21,*) nsteps,inw,inlines        ! Nage Nwave Nwave_line
               read(21,*) (w(i),i=1,inw)            ! wavelengthes
               read(21,*) (wline(i),i=1,inlines)    ! wave lines
	       close(21)
	       if (inw.gt.imw) imw=inw
ccccccccccccccccccccc
c  if GISSEL STOCHASTIC MODELS (BC_STOCH)  
            elseif (typ(1:lnblnk(typ)).eq."BC_STOCH") then
               libphys='YES'
               write(UO,*) 'reading stochastic library:',
     >                       name2(1:lnblnk(name2))
               j=0
              open (21,file=name2,form='unformatted',status='old',err=4)
	       read (21)  inw,(ws(i),i=1,inw)
c               write(*,*) inw,ws(1),ws(1000),ws(inw)
	       if (inw.gt.imw) imw=inw
               do while (j .lt. 25000)
                  read(21)  tform,gamma,zmet,tauv0,mu,nburst,mstr1,
     >                      mstr0,tlastburst,(fburst(i),i=1,5),
     >                      (ftot(i),i=1,5),age_wm,age_wr          ! 1st param list 
                  read(21)    !  nsage,(ages(i),sfrs(i),i=1,nsage),(sfrav(i),i=1,5)     ! 2nd param list 
                  read(21)    ! (fprop(i),fprop0(i),i=1,inw)            ! flux 
c          write(*,*) j,tform,gamma,zmet,tauv0,mu,nburst,mstr1,age_wr
                  if (tform.ge.agemin .AND. tform.le.agemax)
     >              jstoch=jstoch+1 
                  j=j+1
               enddo
               close(21)    
               write(UO,*) '     nmod =',jstoch,' imw max =',imw 
	     endif  
	   endif  
      enddo
 2    close(1) 
      if (libphys(1:1).eq."Y")  nmod=jstoch 
c
ccccccccccccccccccccc
c  nber of wavelengthes must be lower than 8000 elements
	if (imw.gt.8000) then 
	   write(UO,*) 'maximal number of elements in lambda ',imw
	   write(UO,*) ' must be lower than 8000 --> STOP '
	   goto 10
	endif   

      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G")
     >  phys= zpwork(1:lnblnk(zpwork)) // '/lib_bin/' 
     >       // lib(1:lnblnk(lib)) // '.phys'
        doc = zpwork(1:lnblnk(zpwork)) // '/lib_bin/' 
     >       // lib(1:lnblnk(lib)) // '.doc'
        lib = zpwork(1:lnblnk(zpwork)) // '/lib_bin/' 
     >       // lib(1:lnblnk(lib)) // '.bin'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   INFO ON SCREEN AND DOC
        open(3,file=doc,status='unknown')
        write(UO,'(A)') "############################################"
        write(UO,'(A)') "#  It s translating SEDs to binary library #"
        write(UO,'(A)') "#     with the following options :          "
	write(UO,'(2A)')"# Config file     : ",config(1:lnblnk(config))
	write(UO,'(2A)')"# Library type    : ",sedtyp(1:1)
	write(UO,'(A,1x,i6)') "# Number of SEDs  :", nmod
c 
	write(3,'(2A)')"CONFIG_FILE      ",config(1:lnblnk(config))
	write(3,'(2A)')"LIB_TYPE         ",sedtyp(1:1)
	write(3,'(2A)')"LIB_FMT          ",typ(1:lnblnk(typ))
c
	if (sedtyp(1:1).eq."s" .or. sedtyp(1:1).eq."S") then
	   write(UO,'(2A)')     "# STAR_SED    :",mod(1:lnblnk(mod))
	   write(UO,'(2A)')     "# STAR_LIB    :",lib(1:lnblnk(lib))
	   write(UO,'(2A)')     "# STAR_LIB doc:",doc(1:lnblnk(doc))
	   write(UO,'(A,I8)')   "# STAR_WMAX   :",imw
	   write(UO,'(A,E12.6)')"# STAR_FSCALE :",fscale
c
	   write(3,'(2A)')     "STAR_SED    ",mod(1:lnblnk(mod))
	   write(3,'(2A)')     "STAR_LIB    ",lib(1:lnblnk(lib))
	   write(3,'(A,I8)')   "STAR_WMAX   ",imw
	   write(3,'(A,E12.6)')"STAR_FSCALE ",fscale
	elseif (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q") then
	   write(UO,'(2A)')     "# QSO_SED    :",mod(1:lnblnk(mod))
	   write(UO,'(2A)')     "# QSO_LIB    :",lib(1:lnblnk(lib))
	   write(UO,'(2A)')     "# QSO_LIB doc:",doc(1:lnblnk(doc))
	   write(UO,'(A,I8)')   "# QSO_WMAX   :",imw
	   write(UO,'(A,E12.6)')"# QSO_FSCALE :",fscale
c
	   write(3,'(2A)')     "QSO_SED    ",mod(1:lnblnk(mod))
	   write(3,'(2A)')     "QSO_LIB    ",lib(1:lnblnk(lib))
	   write(3,'(A,I8)')   "QSO_WMAX   ",imw
	   write(3,'(A,E12.6)')"QSO_FSCALE ",fscale
	elseif (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
	   write(UO,'(2A)')     "# GAL_SED    :",mod(1:lnblnk(mod))
	   write(UO,'(2A)')     "# GAL_LIB    :",lib(1:lnblnk(lib))
	   write(UO,'(2A)')     "# GAL_LIB doc:",doc(1:lnblnk(doc))
	   write(UO,'(2A)')     "# GAL_LIB phys:",phys(1:lnblnk(phys))
	   write(UO,'(2A)') "# SEL_AGE    :",fileage(1:lnblnk(fileage))
	   write(UO,'(A,I8)')   "# GAL_WMAX   :",imw
	   write(UO,'(A,E12.6)')"# GAL_FSCALE :",fscale
           write(UO,'(A,1x,2(E12.6,1x))') "# AGE_RANGE   ",agemin,agemax
c
	   write(3,'(2A)')     "GAL_SED    ",mod(1:lnblnk(mod))
	   write(3,'(2A)')     "GAL_LIB    ",lib(1:lnblnk(lib))
	   write(3,'(A,I8)')   "GAL_WMAX   ",imw
	   write(3,'(A,E12.6)')"GAL_FSCALE ",fscale
	   write(3,'(2A)')     "SEL_AGE    ",fileage(1:lnblnk(fileage))
           write(3,'(A,1x,2(E12.6,1x))') "AGE_RANGE   ",agemin,agemax
	endif   
        write(UO,'(A)')"#######################################"
c      STOP
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  initialisation
        do i = 1, 500
	  tb(i-1)  =0.
	enddo  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Writing in output binary & doc files
cccccccccccccccccccccccccccccccccccccccccccccccccccc 
	write(UO,*) ' writing in binary library ...'
c   
        open(1,file=mod,status='old',err=41)
c       Ouput: record length per row in bit (*8) with imw: nb max of lambda
        reclmax = 8*(imw+4)
c       Opening output binary file  
        open(4,file=lib,form='unformatted',status='unknown',
     &  access='direct',recl=reclmax)

        if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
c     Opening output binary Phys : max 10 parameters 
          recpmax=8*(10)
c         open(5,file=phys,form='unformatted',status='unknown',
c     &   access='direct',recl=recpmax)
          open(5,file=phys,status='unknown')
        endif
c     ASCII  doc file associated
        write(3,'(A,I6)')  "NUMBER_SED ",nmod
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      if (libphys(1:1) .eq. 'N') then 
      nrec   = 0
      dummy  = -1
      j      = 0
      jstoch = 0 
c 
      do while (.true.)               !  loop on model list 
c   Read sed file and type from mod
         read(1,'(A)',end=12) str 
	 call val_string(str,paravc,test)
	 name=paravc(1)
c  check that SED is used 
         if (str(1:1).eq.'#' .OR. str(1:1).eq.' ' 
     >      .OR. str(1:1).eq.char(9)
     >      .OR. str(1:1).eq.char(13)) goto 11
c  model ok 
	 j = j + 1 
         d4000    = -99.      !   D4000 = Int4050-4250 Fl / Int 3750_3950 dL   dimensionless
         tau      = -99.      !   Tau (yr)   only for SFH as exp-t/Tau
         age      = -99.      !   Age  (yr)
         mass     = -99.      !   mass (Mo)
         sfr      = -99.      !   SFR (Mo/yr)
         zmet     = -99.      !   Metallicity 
         ltir     = -99.      !   Int_8um^1000um    L_lbda . dlbda    in Log unit Lo
         lnir     = -99.      !   Int_2.1um^2.3um L_lbda.dlbda        in Log unit erg/s/Hz 
         luv      = -99.      !   Int_0.21um^0.25um L_lbda . dlbda    in Log unit erg/s/Hz 
         lopt     = -99.      !   Int_0.55um^0.65um L_lbda . dlbda    in Log unit erg/s/Hz
         kmax = 0
         if (test .eq. 1) then 
            typ = 'A'
         else
  	    typ=paravc(2)
         endif
         name2=zpdir(1:lnblnk(zpdir)) // path(1:lnblnk(path))
     >                 // name(1:lnblnk(name))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc     If  ASCII files  (A).......
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 if (typ(1:1) .eq. 'A' ) then 
	      if (sedtyp(1:1).eq."s" .or. sedtyp(1:1).eq."S")
     >                     jtype=j
c     >                     jtype=-1*j
    
	      if (sedtyp(1:1).eq."q" .or. sedtyp(1:1).eq."Q")
     >                     jtype=j
c     >                     jtype=1000+j
	      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G")
     >                     jtype=j
              nsteps   =  1
	      nage     =  1
	      newage(1)= 0.
	      open (2,file=name2,status='old',err=4)
	      i = 0
              do while (.true.)
                i = i + 1
 	        read(2,*,end=32)  w(i),h(i)
		if (h(i).gt.0) then 
		   h(i) = h(i)*fscale
		else
		   h(i) = 0.
		endif
              enddo
 32	      inw = i - 1
              close(2)     
	      write (UO,400)  name,inw,char(13)
	      call flush(UO)
c
c  compute Luminosity  Int(Llbda. dlbda) in UV Opt FIR 
	      if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") then
	         fluxconv = 1./(4*3.141593*100*pc**2)
                 do i = 1, inw 
                    wir(i) = w(i) /10000. 
                    hir(i) = h(i) /fluxconv
                 enddo    
                 if (wir(1).le.8.d0 .and. wir(inw).gt.1000.d0) then
                    call TRAPZD(wir,hir,inw,8d0,1000d0,ltir)
                    if (ltir .gt.0) then 
                       ltir=dlog10(ltir*10000./Lsol)
                    else 
                       ltir=-99.
                    endif
                 endif
                 if (wir(1).le.0.21 .and. wir(inw).gt.0.25) then
                    call TRAPZD(wir,hir,inw,0.21d0,0.25d0,luv)
                    if (luv.gt.0) then 
                         luv=dlog10(luv*10000./400*2300**2/c)
                    else
                         luv=-99.
                    endif
                 endif
                 if (wir(1).le.2.1 .and. wir(inw).gt.2.3) then
                    call TRAPZD(wir,hir,inw,2.1d0,2.3d0,lnir)
                    if (lnir.gt.0) then 
                         lnir=dlog10(lnir*10000./2000.*22000**2/c)
                    else
                         lnir=-99.
                    endif
                 endif
                 if (wir(1).le.0.55 .and. wir(inw).gt.0.65) then
                    call TRAPZD(wir,hir,inw,0.55d0,0.65d0,lopt)
                    if (lopt.gt.0) then 
                         lopt=dlog10(lopt*10000./1000.*6000**2/c)
                    else
                       lopt=-99.
                    endif
                 endif
                 if (wir(1).le.0.375 .and. wir(inw).gt.0.425) then
                    call TRAPZD(wir,hir,inw,0.375d0,0.395d0,ltir0)
                    call TRAPZD(wir,hir,inw,0.405d0,0.425d0,ltir1)
                    if (ltir0.gt.0 .and. ltir1.gt.0) then 
                         d4000=ltir1/ltir0
                    else
                         d4000=-99.
                    endif
                 endif

              endif
c   write lambda then flux in bin output 
              nrec = nrec + 1 
              write(4,rec=nrec) nrec,jtype,dummy,inw,(w(i),i=1,inw)

	      if (sedtyp(1:1).eq."G")  write(5,700) nrec,jtype,dummy,
     >        luv,lopt,lnir,ltir,mass,sfr,zmet,tau,d4000

              nrec = nrec + 1    
              write(4,rec=nrec) nrec,jtype,age,inw,(h(i),i=1,inw)

 	      if (sedtyp(1:1).eq."G")  write(5,700) nrec,jtype,age,
     >        luv,lopt,lnir,ltir,mass,sfr,zmet,tau,d4000

              write(3,'(A,I6.6,3x,i6,2x,i4,1x,i4,1x,A)')
     >            "MOD_",j,jtype,inw,nsteps,name(1:lnblnk(name))
c
ccccccccccccccccccccccccccccccccccccccc
c  If  ASCII/Far IR  files  (LW) : replace Age by Log(Ltir[8-1000um]) 
ccccccccccccccccccccccccccccccccccccccc
         elseif(typ(1:2) .eq.'LW' ) then
              jtype=j
              nsteps    =  1
	      nage      =  1
	      newage(1) =  0.
	      open (2,file=name2,status='old',err=4)
	      i = 0
              do while (.true.)
                i = i + 1
 	        read(2,*,end=33)  w(i),h(i)
		if (h(i).gt.0) then 
		   h(i) = h(i)*fscale
		else
		   h(i) = 0.
		endif
              enddo
 33	      inw = i - 1
              close(2)     
	      write (UO,400)  name,inw,char(13)
	      call flush(UO)
c  compute Luminosity  Int(Llbda. dlbda) in UV Opt FIR 
c   return to Lo/A  by dividing by fluxconv {erg/s/A/cm^2}= 3.197e-07
c             = Lo*fluxconv=Lo/4 PI D^2
              do i = 1, inw 
                   wir(i) = w(i) /10000. 
                   hir(i) = h(i) /(3.197e-7)
              enddo    
              if (wir(1).le.8.d0 .and. wir(inw).gt.1000.d0) then
                  call TRAPZD(wir,hir,inw,8d0,1000d0,ltir)
                  if (ltir .gt.0) then 
                      ltir=dlog10(ltir*10000.)
                  else
                      ltir=-99.
                  endif
              endif
              if (wir(1).le.0.21 .and. wir(inw).gt.0.25) then
                  call TRAPZD(wir,hir,inw,0.21d0,0.25d0,luv)
                  if (luv.gt.0) then 
                      luv=dlog10(luv*10000./400*2300**2/c*Lsol)
                  else
                      luv=-99.
                  endif
              endif
              if (wir(1).le.0.55 .and. wir(inw).gt.0.65) then
                  call TRAPZD(wir,hir,inw,0.55d0,0.65d0,lopt)
                  if (lopt.gt.0) then
                     lopt=dlog10(lopt*10000./1000.*6000**2/c*Lsol)
                  else
                     lopt=-99.
                  endif
              endif
              if (wir(1).le.2.1 .and. wir(inw).gt.2.3) then
                  call TRAPZD(wir,hir,inw,2.1d0,2.3d0,lnir)
                  if (lnir.gt.0) then 
                       lnir=dlog10(lnir*10000./2000.*22000**2/c*Lsol)
                  else
                       lnir=-99.
                  endif
              endif
              if (wir(1).le.0.375 .and. wir(inw).gt.0.425) then
                    call TRAPZD(wir,hir,inw,0.375d0,0.395d0,ltir0)
                    call TRAPZD(wir,hir,inw,0.405d0,0.425d0,ltir1)
                    if (ltir0.gt.0 .and. ltir1.gt.0) then 
                         d4000=ltir1/ltir0
                    else
                         d4000=-99.
                    endif
              endif

c   write lambda then flux in bin output 
              nrec = nrec + 1 
              write(4,rec=nrec) nrec,jtype,dummy,inw,(w(i),i=1,inw)
          write(5,700) nrec,jtype,dummy,luv,lopt,lnir,ltir,
     >                 mass,sfr,zmet,tau,d4000
              nrec = nrec + 1    
              write(4,rec=nrec) nrec,jtype,ltir,inw,(h(i),i=1,inw)
          write(5,700) nrec,jtype,age,luv,lopt,lnir,ltir,
     >                 mass,sfr,zmet,tau,d4000
c
	      write(3,'(A,I6.6,3x,i6,2x,i4,1x,i4,1x,A,2x,F9.6)')
     >            "MOD_",j,jtype,inw,nsteps,name(1:lnblnk(name)),ltir
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  If binary from GISSEL (B) --> GISSEL convention ...
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	 elseif (typ(1:lnblnk(typ)) .eq. "B") then
c  conversion  from Lum {Lo/A} to flux {erg/cm2/s/A}
c  Fluxconv {erg/s/A/cm^2}= Lo{erg/s/A}/[4.pi.(10pc)^2{cm^2}]=3.197e-07 
	      fluxconv = Lsol/(4*3.141593*100*pc**2)
c
	      open (2,file=name2,form='unformatted',status='old',err=4)
	      read  (2)  nsteps,(tbs(i),i=0,nsteps-1)       
	      read  (2)  inw,(ws(i),i=1,inw)
	      do i = 0 , nsteps-1
		 tb(i)=DBLE(tbs(i))
	      enddo
	      do i =1,inw
		 w(i) = DBLE(ws(i))
	      enddo
	      if ( fileage(1:4) .ne. 'NONE') then
                 call CLOSEAGE(fileage,tb,nsteps,cage,nage)
	      else
		 do n = 1, nsteps
		    cage(n-1)=1
		 enddo   
		 nage=nsteps
	      endif 
c    write lambda	 
              nrec = nrec + 1
              write(4,rec=nrec) nrec,j,dummy,inw,(w(i),i=1,inw)
              write(5,700) nrec,j,dummy,luv,lopt,lnir,ltir,
     >                     mass,sfr,zmet,tau,d4000
              k = 0 
 	      do n=1,nsteps
	        read  (2) inw,(hgs(i),i=1,inw)
                if (cage(n-1).eq.1  .and. 
     >              tb(n-1).ge.agemin .and. tb(n-1).le. agemax ) then
                  do i = 1,inw
		     if (hgs(i).lt.0) then
			hg(i) = 0
		     else
			h(i)  = hgs(i)*fluxconv    ! GISSEL conversion in flux (erg/s/cm2/A)
  		        hg(i) = h(i)  * fscale     ! parameter scaling 
		     endif
                  enddo
c        compute Luminosity  Int(Llbda. dlbda) in UV Opt FIR 
                  do i = 1, inw 
                     wir(i) = w(i) /10000. 
                     hir(i) = hg(i)/fluxconv*Lsol  ! erg/s/A
                  enddo    
                  if (wir(1).le.8.d0 .and. wir(inw).gt.1000.d0) then
                     call TRAPZD(wir,hir,inw,8d0,1000d0,ltir)
                     if (ltir .gt.0) then 
			ltir=dlog10(ltir*10000./Lsol)
                     else
			ltir=-99.
		     endif 
                  endif
                  if (wir(1).le.0.21 .and. wir(inw).gt.0.25) then
                    call TRAPZD(wir,hir,inw,0.21d0,0.25d0,luv)
                    if (luv.gt.0) then 
		       luv=dlog10(luv*10000./400*2300**2/c)
                    else
		       luv=-99.
		    endif
                  endif
                  if (wir(1).le.0.55 .and. wir(inw).gt.0.65) then
                    call TRAPZD(wir,hir,inw,0.55d0,0.65d0,lopt)
                    if (lopt.gt.0) then
		       lopt=dlog10(lopt*10000./1000.*6000**2/c)
                    else
		       lopt=-99.
		    endif
                  endif
                  if (wir(1).le.2.1 .and. wir(inw).gt.2.3) then
                    call TRAPZD(wir,hir,inw,2.1d0,2.3d0,lnir)
                    if (lnir.gt.0) then 
                         lnir=dlog10(lnir*10000./2000.*22000**2/c)
                    else
                         lnir=-99.
                    endif
		 endif
                 if (wir(1).le.0.375 .and. wir(inw).gt.0.425) then
                    call TRAPZD(wir,hir,inw,0.375d0,0.395d0,ltir0)
                    call TRAPZD(wir,hir,inw,0.405d0,0.425d0,ltir1)
                    if (ltir0.gt.0 .and. ltir1.gt.0) then 
                         d4000=ltir1/ltir0
                    else
                         d4000=-99.
                    endif
                 endif
c    write flux for selected age
                  nrec = nrec+1
                  write(4,rec=nrec) nrec,j,tb(n-1),inw,
     >                   (hg(i),i=1,inw)
c    write the phys file 
c    nrec, model, Age, lg(Luv), lg(LR), lg(Lk),  lg(Ltir), Mass, SFR     Zmet Tau    D4000
c                 (yr) (---- Log erg/s/Hz ---)     (Lo)    (Mo)  (Mo/yr)       (yr) 
                  write(5,700) nrec,j,tb(n-1),luv,lopt,lnir,ltir,
     >            mass,sfr,zmet,tau,d4000
                  k = k + 1
                  newage(k) = tb(n-1) 
                endif
  	      enddo                   
              nage=k
	      write (UO,410)     name,inw,nage,char(13)
	      call flush(UO)
              close(2)     
	      write(3,'(A,I6.6,3x,i6,2x,i4,1x,i4,1x,A)')
     >        "MOD_",j,j,inw,nage,name(1:lnblnk(name))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc If GISSEL BC03 with physical parameters  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	 elseif (typ(1:lnblnk(typ)) .eq. "BC03") then
c  conversion  from Lum {Lo/A} to flux {erg/cm2/s/A}: fluxconv = 3.197e-07
	      fluxconv = Lsol/(4*3.141593*100*pc**2)
c
	      open (2,file=name2,form='unformatted',status='old',err=4)
	      read  (2)  nsteps,(tbs(i),i=0,nsteps-1),
     > ml,mu0,iseg,(xx(i),lm(i),um(i),baux(i),cn(i),cc(i),i=1,iseg),
     > totm,totn,avs,jo,tau0,id,tau1,tau2,tau3,tau4,id2,id3,iop,stelib
c
c              write(*,*) 
c              write(*,*) 
c              write(*,*) 
              tau = DBLE(tau0)
	      read  (2)  inw,(ws(i),i=1,inw)
	      do i = 0 , nsteps-1
		 tb(i)=DBLE(tbs(i))
	      enddo
	      do i =1,inw
		 w(i) = DBLE(ws(i))
	      enddo
	      if ( fileage(1:4) .ne. 'NONE') then
                 call CLOSEAGE(fileage,tb,nsteps,cage,nage)
	      else
		 do n = 1, nsteps
		    cage(n-1)=1
		 enddo   
	      endif
c  write lambda
              nrec = nrec + 1
              write(4,rec=nrec) nrec,j,dummy,inw,(w(i),i=1,inw)
c              write(5,rec=nrec) nrec,j,dummy,mass,sfr,ltir
            write(5,700) nrec,j,dummy,luv,lopt,lnir,ltir,
     >                   mass,sfr,zmet,tau,d4000
              k = 0 
 	      do n=1,nsteps
	        read  (2) inw,(hgs(i),i=1,inw),ix,(fx(i),i=1,ix)
c                write(*,*) inw,ix,(fx(i),i=1,5)
                if (cage(n-1).eq.1 .and.
     >              tb(n-1).ge.agemin .and. tb(n-1).le. agemax ) then
                  do i = 1,inw
		     if (hgs(i).lt.0) then
			hg(i) = 0
		     else
			h(i)  = hgs(i)*fluxconv  ! GISSEL conversion in flux
  		        hg(i) = h(i)* fscale     ! parameter scaling 
		     endif
                  enddo
c        compute Luminosity  Int(Llbda. dlbda) in UV Opt FIR 
                  do i = 1, inw 
                     wir(i) = w(i) /10000. 
                     hir(i) = hg(i)/fluxconv*Lsol
                  enddo    
                  if (wir(1).le.8.d0 .and. wir(inw).gt.1000.d0) then
                    call TRAPZD(wir,hir,inw,8d0,1000d0,ltir)
                    if (ltir .gt.0) then 
		       ltir=dlog10(ltir*10000./Lsol)
                    else
		       ltir=-99. 
		    endif
                  endif
                  if (wir(1).le.0.21 .and. wir(inw).gt.0.25) then
                    call TRAPZD(wir,hir,inw,0.21d0,0.25d0,luv)
                    if (luv.gt.0) then 
		       luv=dlog10(luv*10000./400*2300**2/c)
                    else
		       luv=-99.
                    endif
	          endif
                  if (wir(1).le.0.55 .and. wir(inw).gt.0.65) then
                    call TRAPZD(wir,hir,inw,0.55d0,0.65d0,lopt)
                    if (lopt.gt.0) then 
		       lopt=dlog10(lopt*10000./1000.*6000**2/c)
                    else
		       lopt=-99.
		    endif
                  endif
                  if (wir(1).le.2.1 .and. wir(inw).gt.2.3) then
                    call TRAPZD(wir,hir,inw,2.1d0,2.3d0,lnir)
                    if (lnir.gt.0) then 
                         lnir=dlog10(lnir*10000./2000.*22000**2/c)
                    else
                         lnir=-99.
                    endif
                  endif
                  if (wir(1).le.0.375 .and. wir(inw).gt.0.425) then
                    call TRAPZD(wir,hir,inw,0.375d0,0.395d0,ltir0)
                    call TRAPZD(wir,hir,inw,0.405d0,0.425d0,ltir1)
                    if (ltir0.gt.0 .and. ltir1.gt.0) then 
                         d4000=ltir1/ltir0
                    else
                         d4000=-99.
                    endif
                  endif

c    write flux for selected age
                  nrec = nrec+1
                  write(4,rec=nrec) nrec,j,tb(n-1),inw,(hg(i),i=1,inw)
                  k = k + 1
                  newage(k) = tb(n-1) 
                  ltirg(k)  = ltir
                  lnirg(k)  = lnir
                  luvg(k)   = luv
                  loptg(k)  = lopt
                  d4000g(k) = d4000
                endif
  	      enddo                   
              kmax=k   
              nage=kmax             
	      write (UO,410)     name,inw,kmax,char(13)
	      call flush(UO)
c              write(6,*) 'nrec=',nrec
c  BC03 with phys parameters : 12 records after the sed's.
	      read (2) nsteps,(bol(i),i=0,nsteps-1)
	      read (2) nsteps,(mstr(i),i=0,nsteps-1)
	      read (2) nsteps,(sf(i),i=0,nsteps-1)
	      read (2) nsteps,(evf(i),i=0,nsteps-1)
	      read (2) nsteps,(snr(i),i=0,nsteps-1)
	      read (2) nsteps,(pnr(i),i=0,nsteps-1)
	      read (2) nsteps,( bh(i),i=0,nsteps-1)
	      read (2) nsteps,( sn(i),i=0,nsteps-1)
	      read (2) nsteps,( wd(i),i=0,nsteps-1)
	      read (2) nsteps,( rm(i),i=0,nsteps-1)
c              write(6,*) (rm(i),i=180,200)
c  write in phys file
              nrec0 = nrec-kmax
              k=0
              do n =1,nsteps
                if (cage(n-1).eq.1 .and.
     >              tb(n-1).ge.agemin .and. tb(n-1).le. agemax ) then
                  nrec0=nrec0+1
                  k=k+1
c            write(5,rec=nrec0) nrec0,j,tb(n-1),mstr(n-1),sf(n-1),ltir
            write(5,700) nrec0,j,tb(n-1),luvg(k),loptg(k),lnirg(k),
     >                   ltirg(k),mstr(n-1),sf(n-1),zmet,tau,d4000g(k)
                endif                   
              enddo
c              write(6,*) 'nrec=',nrec,' nrec0=',nrec0
              close(2)     
c   Write in ASCII file 'doc'
c              write(3,600)   j,inw,nage,name
	      write(3,'(A,I6.6,3x,i6,2x,i4,1x,i4,1x,A)')
     >        "MOD_",j,j,inw,nage,name(1:lnblnk(name))
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc If ASCII file from PEGASE2 (F) --> PEGASE conventions ...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	 elseif (typ(1:1) .eq. "F") then   
c automatic conversion for PEGASE2 from Lum {erg/s/A} to flux {erg/cm2/s/A}
c            Fluxconv {1/cm^2}= 1/[4.pi.(10pc)^2{cm^2}]
	      fluxconv = 1./(4*3.141593*100*pc**2)
c	      write(6,*) fluxconv	      
c             Ajoute la possibilite de selectionner des age      
	      if ( fileage(1:4) .ne. 'NONE') then
               open(10,file=fileage(1:lnblnk(fileage)),
     >                status='unknown')
               i = 0
               do while (.true.)
              	i = i + 1
              	read(10,*,end=34) ageref(i)   ! age in Gyears
              	ageref(i) = ageref(i)*1.e9   ! age in years
               enddo
 34	       amax = i - 1
               close(10)          
              else 
               do i=1,500
              	 ageref(i) =-999
               enddo
               amax=500
              endif
c
c                   write(*,*) j,tau
	      open (2,file=name2,status='old',err=4)
c     Description of the scenario
              read(2,'(a)') str  
              do while (str(1:10).ne.'**********')
                 read(2,'(a)') str 
		 call val_string(str,paravc,test)
                 if (paravc(1)(1:lnblnk(paravc(1))).eq.'Type'
     >        .AND. paravc(4)(1:lnblnk(paravc(4))).eq.'formation:')then
                 if (paravc(5)(1:lnblnk(paravc(5))).eq.'2') then
                   read(2,'(a)') str 
		   call val_string(str,paravc,test)
		   call check_float(paravc(2),str2)
		   read(str2,'(E13.6)') tau
                   tau=tau*1.e6
c                   write(*,*) 'in ..',j,tau
                   read(2,'(a)') str 
                 endif 
                 endif 
              end do
              read(2,*) nsteps,inw,inlines        ! Nage Nwave Nwave_line
              read(2,*) (w(i),i=1,inw)            ! wavelengthes
              read(2,*) (wline(i),i=1,inlines)    ! wave lines
c write wavelengthes 
              nrec = nrec + 1
              write(4,rec=nrec) nrec,j,dummy,inw,(w(i),i=1,inw)
              write(5,700) nrec,j,dummy,luv,lopt,lnir,ltir,
     >                     mass,sfr,zmet,tau,d4000
c	      
	      nage=0
              do n=1,nsteps                       !  blocks of fluxes per age
c  extracte ages and stellar masses
                 read(2,'(A)') str
		 call val_string(str,paravc,test)

		 call check_float(paravc(1),str2)
		 read(str2,'(E13.6)') tb(n-1)
		 tb(n-1) = tb(n-1)*1.e6                ! age in year 
                 call check_float(paravc(3),str2)
  	         read(str2,'(E13.6)') mass             ! mass in star 
                 call check_float(paravc(8),str2)
  	         read(str2,'(E13.6)') zmet             ! metallicity of gas
c
                 read(2,'(A)') str 
		 call val_string(str,paravc,test)

                 call check_float(paravc(1),str2)
  	         read(str2,'(E13.6)') ltir0            ! Lbol (erg/s)
                 call check_float(paravc(3),str2)
  	         read(str2,'(E13.6)') ltir1            ! Ldust/Lbol
                 ltir=-99
                 if (ltir0.gt.0 .and. ltir1.gt.0) 
     >              ltir=dlog10(ltir1*ltir0/Lsol)      ! Log(Ldust/Lo)

                 call check_float(paravc(4),str2)
  	         read(str2,'(E13.6)') sfr              ! 
                 sfr=sfr/1.e6                          ! SFR Mo/yr
                 read(2,*) (h(i),i=1,inw)
                 read(2,*) (hline(i),i=1,inlines)
c select age 
                 if(fileage(1:4) .ne. 'NONE')then
                   found=0
                   do i=1,amax
                     if(tb(n-1).ge.ageref(i)-10000.and.
     >                  tb(n-1).le.ageref(i)+10000) found=1
                   enddo
                 else
                  found=1
                 endif 
c
		 if (tb(n-1).gt.0 .and. found.eq.1 .and.
     >               tb(n-1).ge.agemin .and. tb(n-1).le.agemax ) then
		   nage=nage+1 
		   newage(nage) = tb(n-1)
		   do i = 1,inw
		      if (h(i).le.0) then 
			 hg(i) = 0
		      else 
			 h(i) = h(i)*fluxconv      ! automatic PEGASE2 conv. in flux
		         hg(i) = h(i)*fscale       ! additional parameter scaling 
 	              endif
		   enddo   
c        compute Luminosity  Int(Llbda. dlbda) in UV Opt FIR 
                   do i = 1, inw 
                      wir(i) = w(i) /10000. 
                      hir(i) = hg(i)/fluxconv
                   enddo    
c                   if (wir(1).le.8.d0 .and. wir(inw).gt.1000.d0) then
c                      call TRAPZD(wir,hir,inw,8d0,1000d0,ltir)
c                      if (ltir .gt.0) then 
c			 ltir=dlog10(ltir*10000./Lsol)
c                      else
c			 ltir=-99.
c		      endif
c                   endif
                   if (wir(1).le.0.21 .and. wir(inw).gt.0.25) then
                      call TRAPZD(wir,hir,inw,0.21d0,0.25d0,luv)
                      if (luv.gt.0) then 
			 luv=dlog10(luv*10000./400*2300**2/c)
                      else
			 luv=-99.
		      endif
                   endif
                   if (wir(1).le.0.55 .and. wir(inw).gt.0.65) then
                      call TRAPZD(wir,hir,inw,0.55d0,0.65d0,lopt)
                      if (lopt.gt.0) then
			 lopt=dlog10(lopt*10000./1000.*6000**2/c)
                      else
			 lopt=-99.
		      endif
                   endif
                   if (wir(1).le.2.1 .and. wir(inw).gt.2.30) then
                      call TRAPZD(wir,hir,inw,2.1d0,2.3d0,lnir)
                      if (lnir.gt.0) then
			 lnir=dlog10(lnir*10000./2000.*22000**2/c)
                      else
			 lnir=-99.
		      endif
                   endif
                  if (wir(1).le.0.375 .and. wir(inw).gt.0.425) then
                    call TRAPZD(wir,hir,inw,0.375d0,0.395d0,ltir0)
                    call TRAPZD(wir,hir,inw,0.405d0,0.425d0,ltir1)
                    if (ltir0.gt.0 .and. ltir1.gt.0) then 
                         d4000=ltir1/ltir0
                    else
                         d4000=-99.
                    endif
                  endif

c    write flux for selected age
		   nrec=nrec+1
                   write(4,rec=nrec) nrec,j,tb(n-1),inw,
     >                   (hg(i),i=1,inw)
            write(5,700) nrec,j,tb(n-1),luv,lopt,lnir,ltir, 
     >                 mass,sfr,zmet,tau,d4000
		 endif 
              end do
	      write(UO,410)     name,inw,nage,char(13)
	      call flush(UO)
              close(2)

c   Write in ASCII file 'doc'
c              write(3,600)   j,inw,nage,name
	      write(3,'(A,I6.6,3x,i6,2x,i4,1x,i4,1x,A)')
     >        "MOD_",j,j,inw,nage,name(1:lnblnk(name))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc if libphys='Y'
         elseif (typ(1:lnblnk(typ)) .eq. "BC_STOCH") then 
  	    fluxconv = Lsol/(4*3.141593*100*pc**2)
 	    nage     =  1
            imw      =  0
            zmet=-999.
            open (2,file=name2,form='unformatted',status='old',err=4)
	    read (2)  inw,(ws(i),i=1,inw)
	    if (inw.gt.imw) imw=inw
	    do i =1,inw
	       w(i) = DBLE(ws(i))
	    enddo
c        write lambda  for the first model only  
            if (nrec.eq.0) then 
               nrec = nrec + 1
               write(4,rec=nrec) nrec,jstoch,dummy,inw,(w(i),i=1,inw)
c must have 33 columns in total
              write(5,700)  nrec,jstoch,dummy,luv,lopt,lnir,ltir,
     >         (zmet,i=1,26)
c
   	       write (UO,'(2(i8,1x),e12.6,1x,i4,1x,2(E12.6,1x),a1,$)') 
     >         nrec,jstoch,dummy,inw,w(1),w(inw)
            endif
            j=0
            do while (j .lt. 25000)
              mass     = -99.
              sfr      = -99.
              ltir     = -99.      !   Int_8um^1000um    L_lbda . dlbda    in unit Lo
              luv      = -99.      !   in erg/s/Hz 
              lopt     = -99.      !   in erg/s/Hz 
              lnir     = -99.      !   in erg/s/Hz 
              lfuv     = -99.      !   Int_0.13um^0.17um L_lbda . dlbda    in unit Lo
              lnuv     = -99.      !   Int_0.21um^0.25um L_lbda . dlbda    in unit Lo
              loptr    = -99.      !   Int_0.55um^0.65um L_lbda . dlbda    in unit Lo
              lnirk    = -99.      !   Int_2.1um^2.3um   L_lbda . dlbda    in unit Lo
              afuv     = -99.      !   FUV attenuation -2.5log(L_obs/L_unextinc)
              anuv     = -99.      !   NUV attenuation -2.5log(L_obs/L_unextinc)
              a_g      = -99.      !   g   attenuation -2.5log(L_obs/L_unextinc)
              aopt     = -99.      !   R   attenuation -2.5log(L_obs/L_unextinc)
              a_z      = -99.      !   z   attenuation -2.5log(L_obs/L_unextinc)
              anir     = -99.      !   K   attenuation -2.5log(L_obs/L_unextinc)
              ldustc   = -99.      !   Ldust from cold stars (Lbda>4000) 
              ldusth   = -99.      !   Ldust from hot  stars (Lbda<4000) 
              d4000    = -99.      !   D4000   
c
c      monochromatic Luminosity in erg/s/Hz :  Lnu = Int( Llb dLbda) lbd^2/(c.Dlbda) 
c
              read(2) tform,gamma,zmet,tauv0,mu,nburst,mstr1,
     $                mstr0,tlastburst,(fburst(i),i=1,5),
     $                (ftot(i),i=1,5),age_wm,age_wr                   ! 1st param list 
            read(2) nsage,(ages(i),sfrs(i),i=1,nsage),(sfrav(i),i=1,5)! 2nd param list 
              read(2) (fprop(i),fprop0(i),i=1,inw)                    ! flux 
c 
              if (tform.lt.agemin .OR. tform.gt.agemax) goto 14
              jstoch=jstoch+1
              eta=DBLE(tform)
c              eta=(tform)
              do i = 1,inw
                 if (fprop(i).lt.0) then
                    hg(i) = 0.d0
                    hg0(i)= 0.d0 
                 else
	            hg(i) = DBLE(fprop(i))*fluxconv*fscale         ! GISSEL+scaling flux conversion 
	            hg0(i) = DBLE(fprop0(i))*fluxconv*fscale       ! GISSEL+scaling flux conversion 
	         endif
	      enddo
c
c    write flux for selected age
              nrec = nrec+1
c              write(*,*) nrec,jstoch,eta,inw,hg(1),hg(100),hg(inw)
c
              write(4,rec=nrec) nrec,jstoch,eta,inw,(hg(i),i=1,inw)
	      write (UO,'(2(i8,1x),e12.6,1x,i4,1x,2(E12.6,1x),a1,$)') 
     >              nrec,jstoch,eta,inw,hg(1),hg(inw),char(13)
c	   write (UO,'(i8,2x,i6,2x,i4,a1,$)')   j,nrec,inw,char(13)
	      call flush(UO)
c
cccccccccccccccccccccccccccccc
c        compute Luminosity  Int(Llbda. dlbda) from  UV-Opt to FIR + Attenuation + Ldust(<4000)/Ld(>4000) 
              do i = 1, inw 
                 wir(i)  = w(i)  / 10000. 
                 hir(i)  = hg(i) / fluxconv*Lsol
                 hir0(i) = hg0(i)/ fluxconv*Lsol
              enddo    
cccccccccccccccccccccccccccccc
              if (wir(1).le.0.21 .and. wir(inw).gt.0.25) then
                 call TRAPZD(wir,hir,inw,0.21d0,0.25d0,luv)
                 if (luv.gt.0) luv=dlog10(luv*10000./400*2300**2/c)
              endif
              if (wir(1).le.0.55 .and. wir(inw).gt.0.65) then
                 call TRAPZD(wir,hir,inw,0.55d0,0.65d0,lopt)
                 if (lopt.gt.0) lopt=dlog10(lopt*10000./1000.*6000**2/c)
              endif
              if (wir(1).le.2.1 .and. wir(inw).gt.2.30) then
                call TRAPZD(wir,hir,inw,2.1d0,2.3d0,lnir)
                if (lnir.gt.0) lnir=dlog10(lnir*10000./2000.*22000**2/c)
              endif
ccccccccccccccccccccccccccccccc
c  Ldust 
              call TRAPZD(wir,hir,inw,wir(1),wir(inw),ltir1)
              call TRAPZD(wir,hir0,inw,wir(1),wir(inw),ltir0)
              if ( ltir0.gt.ltir1) 
     >                ltir=dlog10((ltir0-ltir1)*10000./Lsol)
c  Ldust from hot stars
              call TRAPZD(wir,hir,inw,wir(1),0.4d0,ltir1)
              call TRAPZD(wir,hir0,inw,wir(1),0.4d0,ltir0)
              if ( ltir0 .gt. ltir1)
     >           ldusth=dlog10((ltir0-ltir1)*10000./Lsol)
c  Ldust from cold stars
              call TRAPZD(wir,hir,inw,0.4d0,wir(inw),ltir1)
              call TRAPZD(wir,hir0,inw,0.4d0,wir(inw),ltir0)
              if ( ltir0 .gt. ltir1)
     >           ldustc=dlog10((ltir0-ltir1)*10000./Lsol)

c   L FUV  Lbd = 1500A and Dlbda=400A
              if (wir(1).le.0.13 .and. wir(inw).gt.0.17) then
                 call TRAPZD(wir,hir,inw,0.13d0,0.17d0,lfuv)
                 if (lfuv.gt.0)  lfuv=dlog10(lfuv*10000./Lsol)
	      endif
c  L NUV  Lbda= 2300A  and Dlbda=400A
              if (wir(1).le.0.21 .and. wir(inw).gt.0.25) then
                 call TRAPZD(wir,hir,inw,0.21d0,0.25d0,lnuv)
                 if (lnuv.gt.0)  lnuv=dlog10(lnuv*10000./Lsol)
              endif
c  L R   Lbda=6000A  and Dlbda=1000A
              if (wir(1).le.0.55 .and. wir(inw).gt.0.65) then
                 call TRAPZD(wir,hir,inw,0.55d0,0.65d0,loptr)
                 if (loptr.gt.0) loptr=dlog10(loptr*10000./Lsol)
              endif
c  L K  Lbda=22000A  and Dlbda=2000A
              if (wir(1).le.2.1 .and. wir(inw).gt.2.3) then
                 call TRAPZD(wir,hir,inw,2.1d0,2.3d0,lnirk)
                 if (lnirk.gt.0)  lnirk=dlog10(lnirk*10000./Lsol)
              endif
c A_fuv  
              if (wir(1).le.0.13 .and. wir(inw).gt.0.17) then
                 call TRAPZD(wir,hir,inw,0.13d0,0.17d0,ltir1)
                 call TRAPZD(wir,hir0,inw,0.13d0,0.17d0,ltir0)
                 if (ltir0.gt.ltir1) afuv=-2.5*dlog10(ltir1/ltir0)
              endif
c A_nuv  
              if (wir(1).le.0.21 .and. wir(inw).gt.0.25) then
                 call TRAPZD(wir,hir,inw,0.21d0,0.25d0,ltir1)
                 call TRAPZD(wir,hir0,inw,0.21d0,0.25d0,ltir0)
                 if (ltir0.gt.ltir1) anuv=-2.5*dlog10(ltir1/ltir0)
              endif
c A_g  
              if (wir(1).le.0.43 .and. wir(inw).gt.0.47) then
                 call TRAPZD(wir,hir,inw,0.43d0,0.47d0,ltir1)
                 call TRAPZD(wir,hir0,inw,0.43d0,0.47d0,ltir0)
                 if (ltir0.gt.ltir1) a_g=-2.5*dlog10(ltir1/ltir0)
              endif
c A_r 
              if (wir(1).le.0.58 .and. wir(inw).gt.0.62) then
                 call TRAPZD(wir,hir,inw,0.58d0,0.62d0,ltir1)
                 call TRAPZD(wir,hir0,inw,0.58d0,0.62d0,ltir0)
                 if (ltir0.gt.ltir1) aopt=-2.5*dlog10(ltir1/ltir0)
              endif
c A_z 
              if (wir(1).le.0.83 .and. wir(inw).gt.0.87) then
                 call TRAPZD(wir,hir,inw,0.83d0,0.87d0,ltir1)
                 call TRAPZD(wir,hir0,inw,0.83d0,0.87d0,ltir0)
                 if (ltir0.gt.ltir1) a_z=-2.5*dlog10(ltir1/ltir0)
              endif
c A_k  
              if (wir(1).le.2.1 .and. wir(inw).gt.2.3) then
                 call TRAPZD(wir,hir,inw,2.1d0,2.3d0,ltir1)
                 call TRAPZD(wir,hir0,inw,2.1d0,2.3d0,ltir0)
                 if (ltir0.gt.ltir1) anir=-2.5*dlog10(ltir1/ltir0)
              endif
c d4000  Bruzual 83
              if (wir(1).le.0.375 .and. wir(inw).gt.0.425) then
                 call TRAPZD(wir,hir,inw,0.375d0,0.395d0,ltir0)
                 call TRAPZD(wir,hir,inw,0.405d0,0.425d0,ltir1)
                 if (ltir0.gt.0 .and. ltir1.gt.0) then
                      d4000=ltir1/ltir0
                 else
                      d4000=-99.
                 endif
              endif
c
            write(5,700) nrec,jstoch,tform,luv,lopt,lnir,ltir,mstr1,  !  7 first phys-para
     >   (sfrav(i),i=3,4),gamma,zmet,tauv0,mu,age_wm,age_wr,
     >    tlastburst,(fburst(i),i=3,4),DBLE(nburst),mstr0,
     >    lfuv,lnuv,loptr,lnirk,afuv,anuv,a_g,aopt,a_z,anir,ldusth,d4000
c   sfrav / fburst : 1e6 1e7 1e8 1e9 2e9
c
c   Write in ASCII file 'doc'  : just dont write because to big !!
c	     write(3,'(A,I6.6,3x,i6,2x,i4,1x,i4,1x,A)')
c     >           "MOD_",j,j,inw,nage,name2(1:lnblnk(name2))
 14	      j=j+1
           enddo     
 	   close(2)
           write(UO,*) ' selected  models :',jstoch,' over ',j,nrec
	 endif        ! close if  for type of files (A/B/F/...)
c
 11	 continue 
      enddo           ! close loop for number of models in list  
 12   close (1)
      call flush(UO)
      write(UO,*) ' DONE            '
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	cr_date=fdate()
	write(3,'(2A)')      "CREATION_DATE ",cr_date(1:lnblnk(cr_date))
	write(3,'(A,3x,I8)') "NUMBER_ROWS   ",nrec
	write(3,'(A,3x,I8)') "RECORD_LENGTH ",reclmax
	write(3,'(A,3x,I8)') "NUMBER_LBDA   ",inw
	write(3,'(A,3x,I8)') "NUMBER_AGE    ",nage
c	write(3,200)  (newage(i),i=1,nage) 



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Close files
	close (3)
        close(4)
        if (sedtyp(1:1).eq."g" .or. sedtyp(1:1).eq."G") close(5)
        if(UO.eq.30) close(30)

	goto 10
 4	write (UO,*) 'File ',name2(1:lnblnk(name2)),' not found'
 41	write (UO,*) 'File ',mod(1:lnblnk(mod)),' not found'
c
c100	format (4i10)
c200	format (1p5e15.7)
c201	format (e12.6,1x,1p5e15.7)
c300	format (a50)
400	format (a40,1x,i4," points per record ",a1,$)
410	format (a40,1x,i4," points per record ",i4," nb of age",a1,$)
c500     format (1x,i4,' # number of  models')
c600     format (i4,2x,i4,1x,i4,1x,a50)
 700	format (i8,2x,i6,2x,40(E13.6,2x))
10	end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        SUBROUTINE CLOSEAGE(file,agei,nagei,cage,amax)
c
        implicit none
        integer*4 i,j,cage(0:500),nagei,amax,ib
        real*8    agei(0:500),dage,ageref(500)
        character file*4096
c
	ib=0
c   Read the reference age in file
        open(10,file=file,status='unknown')
        i = 0
        do while (.true.)
           i = i + 1
           read(10,*,end=1) ageref(i)   ! age in Gyears
           ageref(i) = ageref(i)*1.e9   ! age in years
        enddo
 1	amax = i - 1
        close(10)          
c
        do j = 1,nagei
          cage(j) = 0
        enddo
        do i = 1,amax
          dage = 9.e10
          do j = 1,nagei
            if (dabs(agei(j-1)-ageref(i)).le.dage) then
              ib = j -1
              dage = dabs(agei(j-1)-ageref(i))
            endif      
            if (j.eq.nagei) cage(ib) = 1      
          enddo
        enddo
c
        return
        end




