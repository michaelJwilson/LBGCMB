c
      PROGRAM  ZPHOT
c   1 option : zphot 
c                    -c zphot.para       : config file 
c                      
c     last modif 09/11/00
c     Measurement of Phot. redshift
c     Prelimary version with no spectra extraction 
c     Last modif : 1 : add context in input catalg specifying the 
c                      filters used for each object
c                  2 : if prior used : choice of the band 
c                  3 : correction of NaN values if all bands are UPPER_LIMITS
c                  4 : Add Distance Modulus for applying Abs-Mag  
c                  5 : Add Option on command line
c                  6 : Read more than 1 library 
c                  7 : Add fast mode
c                  8 : Add Interpolation of best z 
c                  9 : Add Zproba in various z intervalles
c                  10: Add zform per type (only for synthetical libraries)
      implicit none
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  VARIABLES 
      integer*4 chisize,nbf,inlib,wmax,zadapt
      INCLUDE 'out_unit.decl'
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_zchi.decl'
      INCLUDE 'dim_lib.decl'
      INCLUDE 'dim_wave.decl'      
c      parameter (wmax=8000)
      parameter (zadapt=30000)
      integer*4 imag,imax,imin,nf,i,j,k,nobj,nobjm,imagm
      integer*4 nchi,spec,model,lnblnk
      integer*4 cont,gbcont,bp,bused(nbf),bdincl,babs,bdscal
      integer*4 recb(chisize),recb0(chisize),imasb(chisize),recmin(3)
      integer*4 colibmax,chimax,chitrans,pass
      integer*4 nlib,imasmin(3),ndz,nb(100)
      integer*4 nsp,nspmax,nspmaxg,nspmaxs,nspmaxq
      integer*4 irecz0,numlib,iext,buscal(nbf)
      integer*4 nebv,modext(2),iopa(81),nerr,npdz
      integer*4 modlib(inlib),nmeas(5)
      integer*4 paravi(100),test,nbused,liblength,nmod,nmodg
      integer*4 iwpara,iwout,npara(100)
      integer*4 nzf
      integer*4 nbul,nbus
      integer*4 cat_fmt
      real*8    zf(100),zflib(inlib),zfb(chisize),zfmin(3),zfmod
      real*8    fsp(wmax),fgal(5,wmax),wgal(5,wmax),fq(wmax),wq(wmax)
      real*8    wsp(wmax),fst(wmax),wst(wmax),magabsl(2)
      real*8    extb(chisize),extmin(3),exti,mod_dist,mod_dist2,z0
      real*8    mag_abs(3),mag_absb(chisize),dist_mod(chisize)
      real*8    mod_distb(3),mabs(nbf),mabsq(nbf),kapq(nbf),kap(nbf)
      real*8    extic(2,wmax),opal(81,wmax),opat(81,wmax)
      real*8    h0,om0,l0,dz,abs_mag,lbdmaxwr,lbdminwr
      real*8    ageb(chisize),dmb(chisize),zb(chisize)
      real*8    kcorb(nbf,chisize)
      real*8    zp(nbf),mag(nbf),kcor(nbf),z,age,zs,zstep,zmax
      real*8    ab(nbf),sab(nbf),abo(nbf),sabo(nbf),avmagt,avmago
      real*8    zmin(3),dm,agemin(3),dmmin(3),zintb,zbest,ztemp
      real*8    chi2,chimin(3),chi(2,chisize),chibest
      real*8    xp(chisize),yp(chisize)
      real*8    abcor(nbf),flmoy(nbf),flwidth(nbf)
      real*8    maxlz(chisize),dzml,summl,mlarea,pdz(100)
      real*8    zinf,zsup,zmax_gal
      real*8    int_pdz(100),dzpdz(50)
      real*8    extilib(inlib),maglibf(nbf,inlib),klib(nbf,inlib)
      real*8    zlib(inlib),agelib(inlib),maglib(nbf,inlib)
      real*8    funz,lmasi,lmass,fac_err
      real*8    fcorr(nbf)
c baysian 
      real*8    chibay(chisize),zbay,zmed,zbayi,zbays
      real*8    barea
c prior
      real*8    pweight,iab,dzp
      real*8    nzprior,nzprior2,nzprior3,nzprior4,nzpriorVVDS5
c      real*8    color_rf,bp_B,bp_I
      real*8    timy,zform(100),zform2,tused,tzform,tuniv
      real*8    paravr(100),ebv(100),dz_win,min_thres,min_err(100)
      real*8    dchi,z68i,z68s,z90i,z90s,z99i,z99s
      real*8    val1,val2
c      real*8    pb68,pb90,pb99
      real*8    probz(10),zpdzi(10),zpdzs(10)
      real*8    zml68i,zml68s,zml90i,zml90s,zml99i,zml99s
      character*1024 str_inp,str
      character*512 param
      character*512 paravc(100),str_ch
      character*512 zpdir,zpwork,ospec,config,file
      character*512 filters,colib(5),cat,outf,gallib,qsolib,starlib
      character*512 valf(nbf),magtyp,catmag,zfix,cattyp,extlaw
      character*512 outsp,typm,fileop(100)
      character*512 outpara,wpara(100),str_out(100)
      character*512  cr_date,fdate
c  relative to color-sort
      character*512 sel,fastmod,zintp
      integer*4 numcol,f_index(nbf,inlib),reclist(inlib)
      integer*4 numrec,redoing
      real*8    sigcol,fsort(nbf,inlib),magm(nbf)
c  for absolute magnitudes 
      integer*4  fobs(chisize,nbf,nbf),fobs4(chisize,nbf,nbf)
      real*8     minkcol(chisize,nbf,nbf)
      integer*4  method(nbf),nmeth
      integer*4  bapp(nbf),recmin0(3),goodfilter(nbf),index
      real*8     magm0(nbf),minkcolor(nbf)
      integer*4  magabscont(nbf),macont,mbused(nbf,nbf)
c  Ajout pour methode 4
      integer*4  l,m,nbBinZ,bappOp(100)
      real*8     zbmin(100),zbmax(100)
c  Output file for PDZ and ABS filter
      integer*4      nfabs,ufabs,fabs
      integer*4      pdz_fabs(nbf)
c      integer*4      pdz_cont(nbf),pdz_meth(nbf)
      real*8         pdz_mabsz,pdz_mabs(nbf,chisize),pdz_kcor(nbf)
      real*8         pdz_z,pdz_dm,pdz_distm
      character*512  outpdz,pdz_file,pdz_abs(nbf),pdz_mod,pdz_zph
c
ccccccccccAUTO-ADAPT 
c  Ajout pour l'auto-apdapt
      character*512 autoadapt,adapterror
      integer*4 degre,nbshift,cont_ada(zadapt),mod_ada(zadapt)
      integer*4 meth_ada,ngals_ada,adcont,admmin,admmax
      integer*4 realise,iteration,fl_auto,fl1,fl2,iter_best
      real*8    x,x2,x3,auto_thresmin,auto_thresmax
      real*8    ab_ada(zadapt,nbf),sab_ada(zadapt,nbf)
      real*8    magm_ada(zadapt,nbf),zs_ada(zadapt)
      real*8    a0(nbf),a0in(nbf),a1(nbf),a1in(nbf),a2(nbf),a2in(nbf)
      real*8    a3(nbf),a3in(nbf)
      real*8    a0best(nbf),a1best(nbf),a2best(nbf),a3best(nbf)
      real*8    min_errbest(nbf)
      real*8    residu,res_best,chiin,chifit
      real*8    adzmin,adzmax,corr,shift(nbf)
      real*8    maglibini(nbf,inlib),maglibfini(nbf,inlib)
c common for zadapt 
      common /func_int/ k,ngals_ada,fl1,fl2,fl_auto,meth_ada,imagm
      common /func_tab/ ab_ada,sab_ada,magm_ada,zs_ada
      common /func_tab2/ cont_ada,mod_ada
      common /min_err/min_err
ccccccccccAUTO-ADAPT

ccccccccccTRY C17 METHOD
      integer*4 nbC17,C17fl(100)
      real*8    colC17(nbf),errC17(nbf)
ccccccccccTRY C17 METHOD

ccccccccccccccccccccccccccccccccc
c     ADD pour interdire des filtres 
      integer*4 contforb,new_cont
ccccccccccccccccccccccccccccccccc
c
      EXTERNAL      bdincl 
      external      funz,nzprior,nzprior2,nzprior3,nzprior4,timy

c     If you want to give a name to the output screen file
      if(UO.ne.6)then
           open(UO,file='ZPHOT.screen',status='unknown') 
      endif


c
cccccccccccccccccccccccccc
c Environemental  Variable 
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
      param='zphot'
      call get_help(test)
      if (test .eq. 1) call help(param)
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c Initialisation of Input parameter from config file
       param='-c'
      call get_conf(param,config,test)
      if (test.ne.1) call err_option(param,1)
      call get_path(config)
cccccccc PRIMARY  OPTIONS  ccccccccccccccccc
      param='-CAT_IN'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1) cat = paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1) call err_option(param,1)
      call get_path(cat)
c
      param='-INP_TYPE'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1) typm = paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
c
      param='-CAT_MAG'  
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  catmag=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
c
      param='-CAT_FMT'  
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  then 
         cat_fmt=-1  
         if ( paravc(1)(1:lnblnk(paravc(1)))  .eq. 'MEME' .OR. 
     >        paravc(1)(1:lnblnk(paravc(1)))  .eq. 'meme') cat_fmt=0 
         if ( paravc(1)(1:lnblnk(paravc(1)))  .eq. 'MMEE' .OR.
     >        paravc(1)(1:lnblnk(paravc(1)))  .eq. 'mmee') cat_fmt=1 
         if (cat_fmt.ne.0 .and. cat_fmt.ne.1) call err_option(param,1)
      elseif (test.ne.1) then
         call err_option(param,1)
      endif
c
      param='-FILTER_FILE'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  filters=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  call err_option(param,1)
c
      param='-ZPHOTLIB'  
      call getc_option(param,config,100,paravc,test)
      if (test.ge.1 .and. test.le.3) then
         do i = 1,test
           colib(i)=paravc(i)(1:lnblnk(paravc(i)))
         enddo  
         numlib=test
      elseif (test.gt.3) then
         write(UO,*) 'More than 3 librairies used --> STOP'
         stop
      else
         call err_option(param,1)         
      endif
c
      param='-PARA_OUT'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1) outpara=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1) call err_option(param,1)
      call get_path(outpara)
c
      
ccccccccc  SECONDARY OPTIONS  ccccccccccccccccccc
      param='-CAT_TYPE'  
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  cattyp=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  cattyp='SHORT'
c 
      param='-ERR_SCALE'
      call getf_option(param,config,100,paravr,test)
      nerr=test
      if (test.ge.1) then
         do i = 1,nerr
           min_err(i) =   paravr(i)
         enddo
      else
         nerr=1
         min_err(1)=-1
      endif   
c
      param='-ERR_FACTOR'
      call getf_option(param,config,1,paravr,test)
      if (test.eq.1) fac_err=paravr(1)
      if (test.ne.1) fac_err=1.0
c
      param='-BD_SCALE'
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  bdscal=paravi(1)
      if (test.ne.1)  bdscal=0
c
      param='-GLB_CONTEXT'
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  gbcont=paravi(1)
      if (test.ne.1)  gbcont=-1
c
ccccccccccccccccccccccccccccccccc
cADD juste pour pouvoir interdire certain filtres 
      param='-FORB_CONTEXT'
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  contforb=paravi(1)
      if (test.ne.1)  contforb=-1
c
      param='-MASS_SCALE'  
      call getf_option(param,config,2,paravr,test)
      if (test.eq.2)   lmasi=paravr(1)
      if (test.eq.2)   lmass=paravr(2)
      if (test.ne.2)   lmasi=0.
      if (test.ne.2)   lmass=0.
c
      param='-MAG_ABS'
      call getf_option(param,config,2,paravr,test)
      if (test.eq.2) magabsl(1)=paravr(1)
      if (test.eq.2) magabsl(2)=paravr(2) 
      if (test.ne.2)   magabsl(1)=0.
      if (test.ne.2)   magabsl(2)=0.
c
      param='-MAG_REF'
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  babs=paravi(1)
      if (test.ne.1)  babs=0
c
      param='-ZFIX'  
      call getc_option(param,config,1,paravc,test)      
      if (test.eq.1)   zfix=paravc(1)
      if (test.ne.1)   zfix='NO'
c
c      param='-NZ_PRIOR'  
c      call geti_option(param,config,3,paravi,test)
c      if (test.eq.3)  bp=paravi(1)
c      if (test.eq.3)  bp_B=paravi(2)
c      if (test.eq.3)  bp_I=paravi(3)
c      if (test.ne.3)  then
c        bp=0
c        bp_B=0
c        bp_I=0
c      endif 
      param='-NZ_PRIOR'  
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  bp=paravi(1)
      if (test.ne.1)  bp=0
c
      param='-CAT_OUT'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1) outf=paravc(1)(1:lnblnk(paravc(1)))
      if (test.eq.1) call get_path(outf)
      if (test.ne.1) outf='zphot.out'
c
      param='-SPEC_OUT'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  outsp=paravc(1)(1:lnblnk(paravc(1)))
      if (test.ne.1)  outsp='NO'
c
      param='-DZ_WIN'
      call getf_option(param,config,1,paravr,test)
      if (test.eq.1)  dz_win=paravr(1)
      if (test.ne.1)  dz_win=0.25 
c
      param='-ZMAX_GAL'
      call getf_option(param,config,1,paravr,test)
      if (test.eq.1)  zmax_gal=paravr(1)
      if (test.ne.1)  zmax_gal=99.99 
c
      param='-MIN_THRES'
      call getf_option(param,config,1,paravr,test)
      if (test.eq.1)  min_thres=paravr(1)
      if (test.ne.1)  min_thres=0.1 
c
      param='-FAST_MODE'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  fastmod=paravc(1)
      if (test.ne.1)  fastmod='NO'
c
      param='-COL_NUM'
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  numcol=paravi(1)
      if (test.ne.1)  numcol=3
c
      param='-COL_SIGMA'  
      call getf_option(param,config,1,paravr,test)      
      if (test.eq.1)  sigcol=paravr(1)
      if (test.ne.1)  sigcol=3.
c
       param='-COL_SEL'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  sel=paravc(1)
      if (test.ne.1)  sel='OR'
c
      param='-Z_INTERP'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  zintp=paravc(1)
      if (test.ne.1)  zintp='NO'
c
      param='-PROB_INTZ'
      call getf_option(param,config,100,paravr,test)
      npdz=test
      if (test.ge.1) then
         if (MOD(npdz,2).eq.0) then
           do i = 1,npdz
             int_pdz(i) =   paravr(i)
           enddo
         else
           npdz=1
           int_pdz(1) = 0.
         endif  
      else
         npdz=1
         int_pdz(1) = 0.
      endif   
c
      param='-ZFORM_MIN'
      call getf_option(param,config,100,paravr,test)
      nmod=test      
      if (test.ge.1) then
         do i = 1, nmod
            zform(i)=paravr(i)
         enddo   
      else
         nmod=1
         zform(1)=0.
      endif   
c
      param='-MABS_METHOD'
      call geti_option(param,config,100,paravi,test)
      nmeth=test
      if (nmeth.gt.1 .and. nmeth.le.nbf) then
         do i = 1,nmeth 
            method(i)=paravi(i)
         enddo
      elseif  (nmeth.eq.1) then
         do i = 1,nbf 
            method(i)=paravi(1)
         enddo
      else
         do i = 1,nbf 
            method(i)=0
         enddo
      endif
c
      param='-MABS_CONTEXT'
      call geti_option(param,config,100,paravi,test)
      if (test.gt.1 .and. test.eq.nmeth) then
         do i = 1,nmeth 
            magabscont(i)=paravi(i)
         enddo
      elseif (test.eq.1) then 
        do i = 1,nbf 
            magabscont(i)=paravi(1)
         enddo
      else
         do i = 1,nbf 
            magabscont(i)=-1
         enddo
      endif
c     keyword for MABS_METHOD=2
      param='-MABS_REF'
      call geti_option(param,config,100,paravi,test)
      if (test.gt.1 .and. test.eq.nmeth) then
         do i = 1,nmeth 
            bapp(i)=paravi(i)
         enddo
      elseif (test.eq.1) then
         do i = 1, nbf
            bapp(i)=paravi(1)
         enddo
      else
         do i = 1, nbf
            bapp(i)=i
         enddo
      endif
c     Option of the LF used for MABS_METHOD=4
      param='-MABS_FILT'
      call geti_option(param,config,100,paravi,test)
      call getOpt1Tab(nbBinZ,bappOp,100,test,paravi,param)
c
      param='-MABS_ZBIN'
      call getf_option(param,config,100,paravr,test)
      call getOpt2Tab(nbBinZ,zbmin,zbmax,100,test,paravr,param)
c
ccccccccccAUTO-ADAPT
c
      param='-APPLY_SYSSHIFT'
      call getf_option(param,config,100,paravr,test)
      nbshift=test
      if (test.gt.0) then
        do i = 1,nbshift
           shift(i) = paravr(i)
        enddo
      endif
c
      param='-AUTO_ADAPT'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  autoadapt=paravc(1)
      if (test.ne.1)  autoadapt='NO'
c
      param='-ERROR_ADAPT'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  adapterror=paravc(1)
      if (test.ne.1)  adapterror='NO'
c
      param='-ADAPT_BAND'
      call geti_option(param,config,3,paravi,test)
      if (test.eq.3) then
        fl_auto=paravi(1)
        fl1=paravi(2)
        fl2=paravi(3)
      else
        fl_auto=0
        fl1=0
        fl2=0
      endif
      if(fl_auto.le.0)  autoadapt='NO'
c
      param='-ADAPT_LIM'
      call getf_option(param,config,2,paravr,test)
      if (test.eq.2)then
        auto_thresmin=paravr(1)
        auto_thresmax=paravr(2)
      else
        auto_thresmin=18.0d0
        auto_thresmax=21.5d0
      endif 
c
      param='-ADAPT_POLY'
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1 .and.paravi(1).le.4) then 
         degre=paravi(1)
      else
         degre=1
      endif
c
      param='-ADAPT_METH' 
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  meth_ada=paravi(1)
      if (test.ne.1)  meth_ada=1
c
      param='-ADAPT_CONTEXT'
      call geti_option(param,config,1,paravi,test)
      if (test.eq.1)  adcont=paravi(1)
      if (test.ne.1)  adcont=-1
c
      param='-ADAPT_ZBIN'
      call getf_option(param,config,2,paravr,test)
      if (test.eq.2) adzmin=paravr(1)
      if (test.eq.2) adzmax=paravr(2) 
      if (test.ne.2) adzmin=0.001
      if (test.ne.2) adzmax=6.
c
      param='-ADAPT_MODBIN'
      call geti_option(param,config,2,paravi,test)
      if (test.eq.2) admmin=paravi(1)
      if (test.eq.2) admmax=paravi(2) 
      if (test.ne.2) admmin=1
      if (test.ne.2) admmax=1000     
c
      param='-C17_METHOD'
      call geti_option(param,config,100,paravi,test)
      nbC17=test
      if (nbC17.gt.0) then 
        do i = 1,nbC17
           C17fl(i) = paravi(i)
        enddo
      endif
cccccccccccc PDZ OUTPUT cccccccccccccccccccccccccccc
c
      param='-PDZ_OUT'
      call getc_option(param,config,1,paravc,test)
      if (test.eq.1)  outpdz=paravc(1)
      if (test.ne.1)  outpdz='NONE'
      nfabs=0
      if (outpdz(1:4).ne."NONE" .AND. outpdz(1:4).ne."none") then      
c  filter ref  
        param='-PDZ_MABS_FILT'
        call geti_option(param,config,100,paravi,test)
        if (test.ge.1 ) then
           nfabs=test
           do i = 1, nfabs  
              pdz_fabs(i)=paravi(i)
           enddo
        endif
c  opening Mag abs files 
        do i = 1, nfabs           
           write(ospec,'(i2.2)') pdz_fabs(i) 
           pdz_abs(i)=outpdz(1:lnblnk(outpdz))//'.abs'//ospec(1:2)
        enddo
        pdz_file= outpdz(1:lnblnk(outpdz)) // '.pdz'
        pdz_mod = outpdz(1:lnblnk(outpdz)) // '.mod'
        pdz_zph = outpdz(1:lnblnk(outpdz)) // '.zph'
      endif 
cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  INFO PARAMETERS ON SCREEN 
      write(UO,'(A)')  "#######################################"
      write(UO,'(A)')  "# PHOTOMETRIC REDSHIFT with OPTIONS   #"
      write(UO,'(2A)') "# CAT_IN       : ",cat(1:lnblnk(cat))
      write(UO,'(2A)') "# CAT_OUT      : ",outf(1:lnblnk(outf))
      write(UO,'(2A)') "# PARA_OUT     : ",outpara(1:lnblnk(outpara))
      write(UO,'(2A)') "# INP_TYPE     : ",typm(1:lnblnk(typm))
      write(UO,'(A,2x,I4)') "# CAT_FMT[0:MEME 1:MMEE]: ",cat_fmt
      write(UO,'(2A)') "# CAT_MAG      : ",catmag(1:lnblnk(catmag))
      write(UO,'(2A)') "# FILTER_FILE  : ",filters(1:lnblnk(filters))
      write(UO,'(7A)') 
     >"# ZPHOTLIB     : ",(colib(i)(1:lnblnk(colib(i))),' ',i=1,numlib)
      write(UO,'(A,50(f6.3,1x))')
     >                    "# ERR_SCALE    : ",(min_err(i),i=1,nerr)
      write(UO,'(A,f6.3)') "# ERR_FACTOR   : ",fac_err
      write(UO,'(A,I10)')  "# BD_SCALE     : ",bdscal
      write(UO,'(A,I10)')  "# GLB_CONTEXT  : ",gbcont
      write(UO,'(A,f6.3)') "# ZMAX_GAL     : ",zmax_gal
      write(UO,'(A,f6.3)') "# DZ_WIN       : ",dz_win
      write(UO,'(A,f6.3)') "# MIN_THRES    : ",min_thres
      write(UO,'(A,2f8.3)')"# MASS_SCALE   : ",lmasi,lmass 
      write(UO,'(A,2f8.3)')"# MAG_ABS      : ",magabsl(1),magabsl(2)
      write(UO,'(A,I10)')  "# MAG_REF      : ",babs
      write(UO,'(A,100(f6.2,1x))')
     >                    "# ZFORM (model): ",(zform(i),i=1,nmod)
      write(UO,'(A,I10)')  "# NZ_PRIOR     : ",bp 
      write(UO,'(2A)')     "# Z_INTERP     : ",zintp(1:lnblnk(zintp))
      write(UO,'(A,50(f6.3,1x))')
     >                     "# PROB_INTZ    : ",(int_pdz(i),i=1,npdz)

      if (nmeth.eq.1) then 
        write(UO,'(A,I4)')   "# MABS_METHOD  : ",method(1)
        write(UO,'(A,I10)')   "# MABS_CONTEXT : ",magabscont(1)
        if (method(1).eq.2)
     >  write(UO,'(A,I4)')   "# MABS_REF     : ",bapp(1)
      else
        write(UO,'(A,100(I4,1x))') 
     >   "# MABS_METHOD  : ",(method(i),i=1,nmeth)
        write(UO,'(A,100(I10,1x))') 
     >   "# MABS_CONTEXT : ",(magabscont(i),i=1,nmeth)    
        write(UO,'(A,100(I4,1x))') 
     >   "# MABS_REF     : ",(bapp(i),i=1,nmeth)    
      endif
      test=0
      do k=1,nmeth
        if (method(k).eq.4 .and.test.eq.0) then
          test=1 
          write(UO,'(A,50(I4,1x))')   "# MABS_FILT    : ",
     >  (bappOp(i),i=1,nbBinZ)
          write(UO,'(A,50(f6.3,1x))')"# MABS_ZBIN    : ",
     >  (zbmin(i),zbmax(i),i=1,nbBinZ)
        endif
      enddo

      if (nbshift.gt.0) write(UO,'(A,100(f6.2,1x))')
     > "# APPLY_SYSSHIFT: ",(shift(i),i=1,nbshift)
      write(UO,'(2A)')     "# AUTO_ADAPT   : ",
     >                      autoadapt(1:lnblnk(autoadapt))
      if (autoadapt(1:1) .eq. "Y" .or. autoadapt(1:1).eq."y") then
       write(UO,'(2A)')     "# ERROR_ADAPT  : ",
     >                      adapterror(1:lnblnk(adapterror))
       write(UO,'(A,3(I10,1x))')
     >                    "# ADAPT_BAND   : ",fl_auto,fl1,fl2
       write(UO,'(A,2(f6.2,1x))')"# ADAPT_LIM    : ",auto_thresmin,
     > auto_thresmax
       write(UO,'(A,I10)') "# ADAPT_POLY   : ",degre
       write(UO,'(A,I10)') "# ADAPT_METH   : ",meth_ada
       write(UO,'(A,I10)') "# ADAPT_CONTEXT: ",adcont
       write(UO,'(A,2(f6.3,1x))') "# ADAPT_ZBIN   : ",adzmin,adzmax
       write(UO,'(A,2(I8,1x))')   "# ADAPT_MODBIN : ",admmin,admmax
      endif
      if (nbC17.gt.0) write(UO,'(A,100(i4,1x))')     
     > "# C17_METHOD: ",(C17fl(i),i=1,nbC17)
c
      write(UO,'(2A)') "# ZFIX         : ",zfix(1:lnblnk(zfix))
      write(UO,'(2A)') "# SPEC_OUT     : ",outsp(1:lnblnk(outsp))
c
      write(UO,'(2A)') "# PDZ_OUT      : ",outpdz(1:lnblnk(outpdz))
      if (outpdz(1:4).ne."NONE" .and. outpdz(1:4).ne."none") then
        write(UO,'(2A)')          "# PDZ_file     : ",
     >                   pdz_file(1:lnblnk(pdz_file))
        write(UO,'(2A)') "# PDZ_ZPH     : ",pdz_zph(1:lnblnk(pdz_zph))
        write(UO,'(2A)') "# PDZ_MOD     : ",pdz_mod(1:lnblnk(pdz_mod))
        write(UO,'(A,100(A,1x))')    "# PDZ_MABS     : ",
     >                   (pdz_abs(i)(1:lnblnk(pdz_abs(i))),i=1,nfabs)
        write(UO,'(A,100(1x,I4))')"# PDZ_MABS_FILT: ",
     >                   (pdz_fabs(i),i=1,nfabs)
c        write(UO,'(A,100(1x,I4))')"# PDZ_MABS_CONT: ",
c     >                   (pdz_cont(i),i=1,nfabs)
c        write(UO,'(A,100(1x,I4))')"# PDZ_MABS_METH: ",
c     >                   (pdz_meth(i),i=1,nfabs)

      endif       
c
      write(UO,'(2A)') "# FAST_MODE    : ",fastmod(1:lnblnk(fastmod))
      if (fastmod(1:1).eq."Y" .or. fastmod(1:1).eq."y") then
       write(UO,'(A,I10)') "# COL_NUM      : ",numcol
       write(UO,'(A,f6.2)')"# COL_SIGMA    : ",sigcol
       write(UO,'(2A)')    "# COL_SEL      : ",sel(1:lnblnk(sel))
      endif   



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c READING INPUT LIBRARIES reclmax
      gallib= ' '
      qsolib= ' '
      starlib=' '
      write(UO,'(A,a1,$)') "reading librairies ...",char(13)
      call flush(UO)
      call read_lib(colib,numlib,inlib,zstep,zmax,dz,
     >colibmax,modlib,extilib,zlib,agelib,maglib,klib,
     >nzf,zf,zflib,
     >starlib,qsolib,gallib,extlaw,nebv,ebv,modext,
     >h0,om0,l0,magtyp,valf,imagm,nmodg,fobs,minkcol)
      write(UO,*) 'number of record in libraries: ',colibmax
c
      if (lnblnk(catmag).ne.0 .and. catmag(1:5).ne.magtyp(1:5)) then
         write(UO,*) '  Magnitude types between libraries '
         WRITE(UO,*) '  and INPUT catalog are different   '
         stop
      endif    
c
      write(UO,'(3A,1x,I9,A,I9)') 'number of models in ',
     > gallib(1:lnblnk(gallib)),' : ',nmodg,' & filters : ',imagm
c
c 
c  read cat if type == LONG for length of string_input
      if (cattyp(1:4).eq.'LONG') then
        open(1,file=cat,status='unknown')      
        do while (.true.)                 ! --> closed at the END
           read(1,'(A)',end=18) str 
           call val_string(str,paravc,test)
           if (paravc(1)(1:1) .ne. '#') then  
              j = 2*imagm+3
              if (test.gt.j) then
                cont= test-j
              else
                cont=0
              endif
              goto 18 
          endif
        enddo   
 18     close(1)
c           write(*,*) test,j,cont

      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  reading output parameter
      open(1,file=outpara,status='unknown')
      iwpara = 0 
      do while (.true.)
         read(1,'(A)',end=13) str 
         if (str(1:1) .ne. '#' .AND. str(1:1) .ne. ' ' 
     >  .AND. str(1:1) .ne. char(9) .AND. str(1:1) .ne. char(13)) then
           call val_string(str,paravc,test)
           iwpara = iwpara + 1
           wpara(iwpara) = paravc(1)
c
          if (paravc(1)(lnblnk(paravc(1)):lnblnk(paravc(1))).eq.')')then
            if (iwpara.gt.1) then 
               npara(iwpara)=imagm+npara(iwpara-1)
            else
               npara(iwpara)=imagm
            endif
          elseif (paravc(1)(1:lnblnk(paravc(1))).eq.'STRING_INPUT') then
            if (iwpara.gt.1) then 
               npara(iwpara)=cont+npara(iwpara-1)
            else
               npara(iwpara)=cont
            endif
          else
            if (iwpara.gt.1) then 
               npara(iwpara)=1+npara(iwpara-1)
            else
               npara(iwpara)=1
            endif
          endif
c          
         endif
      enddo
 13   close(1)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(UO,'(A)') "# Output format   "
      val1 = DINT(DBLE(iwpara)/5.)
      if (MOD(iwpara,5).eq.0) then
          k = IDINT(val1) 
      else
          k = IDINT(val1) + 1
      endif
      do i = 1,k
         imin = (i-1)*5 + 1
         imax = imin + 4
         if (imax.gt.iwpara) imax = iwpara
c        write(UO,'(A,1x,4(A," , "),A)') "# ",
c     >  (wpara(j)(1:lnblnk(wpara(j))),j=imin,imax)
      enddo
c
      write(UO,'(A)') "#######################################"
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  COMPUTING THE FILTER CHARACTERISTICS
      write(UO,'(A,a1,$)') "filter characteristics ...",char(13)
      call flush(UO)
      call zeropoint(filters,zp,abcor,flmoy,flwidth,fcorr,imag)
c
      if (imag.ne.imagm) then
         write(UO,*) 'WARNING : Number of filters is not equal between'
         write(UO,*) filters(1:lnblnk(filters)),' and '
     >             ,colib(1)(1:(lnblnk(colib(1))))
         write(UO,*) ' check in ',config(1:lnblnk(config)),
     >                 ' the file name in FILTER_FILE '
         stop
      endif
c scaling errors check
      if (nerr.ne.imag) then
         do i = 1,imag
           min_err(i) = -1.
         enddo
      endif
c   
ccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check zform_min keyword size with number of models
      if (nmod.ne.nmodg) then
c         nmod=nmodg
         do i = 1, nmodg
            zform(i) = 0
         enddo
c         write(UO,'(A)') "zform_min not set properly -> canceled"
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   READ EXTINCTION AND OPACITY  if Output spectra = Y 
      if (outsp(1:1).eq.'Y'.or.outsp(1:1).eq.'y') then 
        if (extlaw(1:4).ne.'NONE' .and. extlaw(1:4).ne.'none') then 
c  read extinction law 
          file = zpdir(1:lnblnk(zpdir)) //'/ext/' 
     >     // extlaw(1:lnblnk(extlaw))
          open(1,file=file,status='unknown')
          i=0
          do while (.true.)
            i=i+1
            read(1,*,end=11) extic(1,i),extic(2,i)
          enddo
 11       iext=i-1
          close(1)
        else
           iext=1
           extic(1,1)=0.
           extic(2,1)=0.
        endif  
c  reading file with  extragalacitic opacity
        file = zpdir(1:lnblnk(zpdir)) // '/opa/OPACITY.dat'
        open(1,file=file,status='unknown')
        do i = 1,81
           read(1,'(A)') str 
           call val_string(str,paravc,test)
           fileop(i)=paravc(2)
        enddo
        close(1)
        do i = 1, 81
           file = zpdir(1:lnblnk(zpdir)) // '/opa/' 
     > // fileop(i)
           open(1,file=file,status='unknown')
           k = 0
           do while (.true.)
              k = k + 1
              read(1,*,end=12) opal(i,k), opat(i,k)
           enddo  
 12        iopa(i) = k - 1
           close(1)
        enddo             
      endif
cccccccccccccccccccccccCHANGEcccccccccccccccccccccccccc
c     Apply systematic shift
      if(nbshift.eq.imagm)then
       do i = 1, colibmax
         do k = 1, imagm
            maglib(k,i)=maglib(k,i)+shift(k)
         enddo
       enddo
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   CONVERSION from theoretical magnitudes to fluxes
      write(UO,'(A,a1,$)') "converting mag to flux ...",char(13)
      call flush(UO)
      do i = 1, colibmax
         do k = 1, imagm
            if (maglib(k,i) .gt. 99) then 
                maglibf(k,i) = 0.
            else
             if (magtyp(1:1).eq.'A') 
     >              maglibf(k,i)=10**(-0.4*(maglib(k,i)+48.60))
             if (magtyp(1:1).eq.'V')
     >              maglibf(k,i)=10**(-0.4*(maglib(k,i)-zp(k)))
            endif
ccccccccccAUTO-ADAPT
            maglibfini(k,i)=maglibf(k,i)
            maglibini(k,i) =maglib(k,i)
ccccccccccAUTO-ADAPT
         enddo          
      enddo   
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c   Sort according to colors 
      if (fastmod(1:1).eq."Y" .or. fastmod(1:1).eq."y") then
        write(UO,'(A)') "sorting the colors..."
        call flush(UO)
        call sort_col(maglibf,imagm,colibmax,numcol,fsort,f_index)
      endif  
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Header for Output file 
c opening files 
      open(4,file=outf(1:lnblnk(outf)),status='unknown')
      if (outpdz(1:4).ne."NONE" .and. outpdz(1:4).ne."none") then
        open(43,file=pdz_zph(1:lnblnk(pdz_zph)),status='unknown') 
        open(44,file=pdz_file(1:lnblnk(pdz_file)),status='unknown')  
        open(45,file=pdz_mod(1:lnblnk(pdz_mod)),status='unknown')  
        do i = 1,nfabs
          ufabs=45+i 
          open(ufabs,file=pdz_abs(i)(1:lnblnk(pdz_abs(i))),
     >          status='unknown')  
        enddo
      endif
      write(4,'(A)') "#######################################"
      write(4,'(A)') "# PHOTOMETRIC REDSHIFT with OPTIONS   #"
      write(4,'(A)') "#          Input catalog              #" 
      write(4,'(2A)') "# CAT_IN      : ",cat(1:lnblnk(cat))
      write(4,'(2A)') "# INP_TYPE    : ",typm(1:lnblnk(typm))
      write(4,'(2A)') "# CAT_MAG     : ",catmag(1:lnblnk(catmag))
      write(4,'(A,2x,I4)') "# CAT_FMT[0:MEME 1:MMEE]: ",cat_fmt
      write(4,'(A)') "#          Input library              #"    
      write(4,'(2A)') "# FILTER_FILE : ",filters(1:lnblnk(filters))
      write(4,'(100A)') "# ",magtyp(1:5),(valf(j)(1:lnblnk(valf(j)))
     >             ,' ',j=1,imagm)
      write(4,'(A,1x,100(f8.3,1x))') "# ",(abcor(j),j=1,imagm)
      write(4,'(8(A))') "# ZPHOTLIB    : ",
     > (colib(i)(1:lnblnk(colib(i))),' ',i=1,numlib)
      write(4,'(8(A))') "# using       : ",gallib(1:lnblnk(gallib))
     >,' ',qsolib(1:lnblnk(qsolib)),' ',starlib(1:lnblnk(starlib))
      write(4,'(2A)') "# FAST_MODE   : ",fastmod(1:lnblnk(fastmod))
      if (fastmod(1:1).eq."Y" .or. fastmod(1:1).eq."y") then
         write(4,'(A,I10)') "# COL_NUM     : ",numcol
         write(4,'(A,f6.2)') "# COL_SIGMA   : ",sigcol
         write(4,'(2A)') "# COL_SEL     : ",sel(1:lnblnk(sel))
      endif
      write(4,'(A,1x,3f6.2)') "# Z_STEP      : ",zstep,zmax,dz
      write(4,'(A,f6.3)')     "# ZMAX_GAL    : ",zmax_gal
      write(4,'(2A)') "# Z_INTERP     : ",zintp(1:lnblnk(zintp))
      write(4,'(A,1x,3f6.2)') "# COSMOLOGY   : ",h0,om0,l0
      write(4,'(2A)') "# EXTINC_LAW  : ",extlaw(1:lnblnk(extlaw))
      write(4,'(A,1x,10f6.3)') "# EB_V : ",(ebv(i),i=1,nebv)
      write(4,'(A,1x,2I8)') "# MOD_EXTINC  : ",(modext(i),i=1,2)
      write(4,'(A,100(f6.3,1x))') "# ERR_SCALE   : ",(min_err(i),
     >  i=1,nerr)
      write(4,'(A,f6.3)') "# ERR_FACTOR  : ",fac_err
      write(4,'(A)') "#          Options                    #"
      write(4,'(A,I10)')  "# BD_SCALE    : ",bdscal
      write(4,'(A,I10)')  "# GLB_CONTEXT : ",gbcont
      write(4,'(A,f6.3)') "# DZ_WIN      : ",dz_win
      write(4,'(A,f6.3)') "# MIN_THRES   : ",min_thres
      write(4,'(A,2(f6.2,1x))') "# MASS_SCALE  : ",lmasi,lmass 
      write(4,'(A,2(f8.2,1x))') "# MAG_ABS     : ",(magabsl(i),i=1,2)
      write(4,'(A,1x,I8)')      "# MAG_REF     : ",babs
      write(4,'(A,100(f6.3,1x))')
     >                    "# PROB_INTZ    : ",(int_pdz(i),i=1,npdz)
      if (nmod .eq. nmodg) then
         write(4,'(A,100(f6.2,1x))') "# ZFORM_MIN(model): ",
     > (zform(i),i=1,nmod)
      endif   

      if (nmeth.eq.1) then 
        write(4,'(A,I4)')   "# MABS_METHOD  : ",method(1)
        write(4,'(A,I10)')   "# MABS_CONTEXT : ",magabscont(1)
        if (method(1).eq.2)
     >  write(4,'(A,I4)')   "# MABS_REF     : ",bapp(1)
      else
        write(4,'(A,100(I4,1x))') 
     >   "# MABS_METHOD  : ",(method(i),i=1,nmeth)
        write(4,'(A,100(I10,1x))') 
     >   "# MABS_CONTEXT : ",(magabscont(i),i=1,nmeth)    
        write(4,'(A,100(I4,1x))') 
     >   "# MABS_REF     : ",(bapp(i),i=1,nmeth)    
      endif
      test=0
      do k=1,nmeth
        if (method(k).eq.4 .and.test.eq.0) then
          test=1 
          write(4,'(A,100(I4,1x))')   "# MABS_FILT    : ",
     >  (bappOp(i),i=1,nbBinZ)
          write(4,'(A,100(f6.3,1x))')"# MABS_ZBIN    : ",
     >  (zbmin(i),zbmax(i),i=1,nbBinZ)
        endif
      enddo

      if (nbshift.gt.0) write(4,'(A,100(f6.2,1x))')
     > "# APPLY_SYSSHIFT: ",(shift(i),i=1,nbshift)
      write(4,'(2A)')     "# AUTO_ADAPT   : ",
     >                      autoadapt(1:lnblnk(autoadapt))
      if (autoadapt(1:1) .eq. "Y" .or. autoadapt(1:1).eq."y") then
       write(4,'(2A)')     "# ERROR_ADAPT  : ",
     >                      adapterror(1:lnblnk(adapterror))
       write(4,'(A,3(I10,1x))')
     >                    "# ADAPT_BAND   : ",fl_auto,fl1,fl2
       write(4,'(A,2(f6.2,1x))')"# ADAPT_LIM    : ",auto_thresmin,
     > auto_thresmax
       write(4,'(A,I10)') "# ADAPT_POLY   : ",degre
       write(4,'(A,I10)') "# ADAPT_METH   : ",meth_ada
       write(4,'(A,I10)') "# ADAPT_CONTEXT: ",adcont
       write(4,'(A,2(f6.3,1x))') "# ADAPT_ZBIN   : ",adzmin,adzmax
       write(4,'(A,2(I8,1x))')   "# ADAPT_MODBIN : ",admmin,admmax
      endif
      if (nbC17.gt.0) write(4,'(A,100(i4,1x))')
     > "# C17_METHOD: ",(C17fl(i),i=1,nbC17)
      write(4,'(A,1x,I8)') "# NZ_PRIOR    : ",bp 
      write(4,'(2A)')      "# ZFIX        : ",zfix(1:lnblnk(zfix))
      write(4,'(2A)')      "# SPEC_OUT    : ",outsp(1:lnblnk(outsp))
c
      write(4,'(2A)') "# PDZ_OUT      : ",outpdz(1:lnblnk(outpdz))
      if (outpdz(1:4).ne."NONE" .and. outpdz(1:4).ne."none") then
        write(4,'(2A)')          "# PDZ_file     : ",
     >                   pdz_file(1:lnblnk(pdz_file))
        write(4,'(2A)') "# PDZ_ZPH     : ",pdz_zph(1:lnblnk(pdz_zph))
        write(4,'(2A)') "# PDZ_MOD     : ",pdz_mod(1:lnblnk(pdz_mod))
        write(4,'(A,100(A,1x))')    "# PDZ_MABS     : ",
     >                   (pdz_abs(i)(1:lnblnk(pdz_abs(i))),i=1,nfabs)
        write(4,'(A,100(1x,I4))')"# PDZ_MABS_FILT: ",
     >                   (pdz_fabs(i),i=1,nfabs)
      endif       
c
      write(4,'(2A)')      "# CAT_OUT     : ",outf(1:lnblnk(outf))
      write(4,'(2A)')      "# PARA_OUT    : ",outpara(1:lnblnk(outpara))
      cr_date=fdate()
      write(4,'(2A)') "# CREATION_DATE: ",cr_date(1:lnblnk(cr_date)) 
      write(4,'(A)')  "# Output format                       #"
      val1 = DINT(DBLE(iwpara)/5.)
      if (MOD(iwpara,5).eq.0) then
          k = IDINT(val1)
      else
          k = IDINT(val1) + 1
      endif
      do i = 1,k
         imin = (i-1)*5 + 1
         imax = imin + 4
         if (imax.gt.iwpara) imax = iwpara
        write(4,'(A,1x,5(A,1x,i3," , "))') "# ",
     >  (wpara(j)(1:lnblnk(wpara(j))),npara(j),j=imin,imax)

      enddo
      write(4,'(A)') "#######################################"
c
c  convert mass in dexp if used
      if (lmasi.gt.0 .and. lmass.gt.0 ) then
         lmasi=10**(lmasi)
         lmass=10**(lmass)
      endif   
c
c Initialisation of chi2 evolution with z and chimin
      if (zmax.le.6.0)  then 
        chimax   = idnint(zmax/zstep) + 1
        chitrans = chimax
      else
        chitrans = idnint(6./zstep) + 1 
        chimax   = chitrans  + idnint((zmax-6.)/dz)
      endif
c Comput the distance modulus for mag_abs
c 
      do k = 1,chimax
        if (k.le.chitrans) then
           chi(1,k) = zstep*(k-1)
        else 
           chi(1,k) = zstep*(chitrans-1) + dz*(k-chitrans)
        endif
        z=chi(1,k)
        dist_mod(k)=funz(z,h0,om0,l0)
      enddo 
c save the z bin for PDZ_OUT 
      if (outpdz(1:4).ne."NONE" .and. outpdz(1:4).ne."none") then
        write(43,'(500(f8.3,1x))') (chi(1,k),k=1,chimax)
      endif
c Number of measured objects
      do k = 1,3 
        nmeas(k)=0
      enddo  
c Initialisation for proba
      do i = 1,3
         zpdzi(i) = 0.
         zpdzs(i) = 0.
         probz(i) = 0.
      enddo   
      zml68i=zpdzi(1)
      zml68s=zpdzs(1)
      zml90i=zpdzi(2)
      zml90s=zpdzs(2)
      zml99i=zpdzi(3)
      zml99s=zpdzs(3)      
c
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     For the case of method 4 in absolute magnitude
c     The imposed observed filter depends on redshift
      test=0
      do j = 1,nmeth
        if(method(j).eq.4 .and. test.eq.0) then
          test=1
          do i = 1, chimax
           do k=1,nbBinZ  !loop en redshift bin
            if(chi(1,i).ge.zbmin(k).and.chi(1,i).lt.zbmax(k)) then
c         If z include in redshift bin k, use the correspondant k filter
             do l=1,imagm
              do m=1,imagm
                fobs4(i,l,m)=bappOp(k)
              enddo
             enddo 
            endif
           enddo
          enddo
        endif
      enddo
c
ccccccccccAUTO-ADAPT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Initialisation of auto adapt
      do k=1,imagm
       a0(k)=0.d0
       a1(k)=0.d0
       a2(k)=0.d0
       a3(k)=0.d0
      enddo
      x=0.d0
      x2=0.d0
      x3=0.d0
      model=0
      chiin=10.d10 
      iteration=0 
      realise=0
      res_best=20000.d0
      if (autoadapt(1:1).eq.'Y'.or. autoadapt(1:1).eq.'y') then
        write(UO,'(A)') "##############################################" 
         write(UO,'(A,1x,I4)') " --> Starting AUTO-ADAPT with method:",
     >   meth_ada
c         write(UO,'(A,a1,$)') " opening file for AUTO-ADAPT ...",char(13)
c         call flush(UO)   
        open(41,file='minuit.dat',status='unknown')
        open(42,file=outf(1:lnblnk(outf))//'.corr',status='unknown')
c        open(42,file='correction.dat',status='unknown')
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     redone from this point if auto-adapt didn't converge
 6    continue
      if (autoadapt(1:1).eq.'Y'.or. autoadapt(1:1).eq.'y') then

c        Reinitialize the number of objects
         do k = 1,3 
           nmeas(k)=0
         enddo  

c        increment iteration
         iteration=iteration+1
         if (realise.lt.1) 
     >      write(UO,'(A,I4)') " --> Iteration :",iteration

c        Stop the iteration if more than 10 and take the values at the best point
         if(iteration.ge.10.and.realise.lt.2) then
            write(UO,*)"MORE THAN 10 ITERATIONS !!"
            write(UO,*)"STOP TRAINING !!"
            realise=2
c           take the best value of iteration
            write(UO,*)"Take the value of iteration :",iter_best
             do k=1,imagm
               a0(k)=a0best(k)
               a1(k)=a1best(k)
               a2(k)=a2best(k)
               a3(k)=a3best(k)
               min_err(k)=min_errbest(k)
            enddo
         endif

c
c        loop on librairy and apply the correction determined with auto-adapt
         do i = 1, colibmax

c         Determine the correction to apply according to the method 
          x=0.d0  
          if (meth_ada.eq.1) then 
c         Predicted color
            if(dabs(maglibini(fl1,i)).le.80.and.
     .         dabs(maglibini(fl2,i)).le.80)then
              x=maglibini(fl1,i)-maglibini(fl2,i)
            endif
          elseif (meth_ada.eq.2) then 
c           Redshift
            x=zlib(i) 
          elseif (meth_ada.eq.3) then          
c            Modele
             x=dble(modlib(i))  
          endif

c         Apply correction to the librairie
          if(degre.ge.3)x2=x**2.d0
          if(degre.ge.4)x3=x**3.d0
          do k=1,imagm 
           corr=a0(k)+a1(k)*x+a2(k)*x2+a3(k)*x3
           maglib(k,i) =maglibini(k,i)+corr
           maglibf(k,i)=maglibfini(k,i)*10.d0**(-0.4*corr)
          enddo

         enddo   
      endif
ccccccccccAUTO-ADAPT


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  OPEN THE PHOTOMETRIC CATALOG  + OUTPUT FILE 
c
c      if(autoadapt(1:1).eq.'N' .or. autoadapt(1:1).eq.'n'
c     >    .or. realise.ge.2 )then
c          write(UO,'(A,a1,$)') 'here we go ... ',char(13)
c         call flush(UO)
c      endif
      irecz0=1
      pweight=1
      cont=  0 
      zs=   -99
      str_inp = ' '
      nobj = 0
      ngals_ada=1
      open(90,file=cat,status='unknown')      
      do while (.true.)                 ! --> closed at the END
         read(90,'(A)',end=8) str 
         call val_string(str,paravc,test)
         if (paravc(1)(1:1) .eq. '#') goto 1 
         read(paravc(1),'(i10)') spec
         do k = 1,imagm
            if (cat_fmt .eq. 0) j = 2*k
            if (cat_fmt .eq. 1) j = k+1
            call check_float(paravc(j),str_ch)
            read(str_ch,'(E24.12)') ab(k)
            if (cat_fmt .eq. 0) j=j+1
            if (cat_fmt .eq. 1) j=j+imagm
            call check_float(paravc(j),str_ch)
            read(str_ch,'(E24.12)') sab(k)
c  keep original values in mag and errmag
            if (typm.eq.'F' .or. typm.eq.'f') then
              if (sab(k).lt.0)  sabo(k) = sab(k) 
              if (sab(k).gt.0) then
                if (ab(k).gt.0) then 
                   sabo(k)=1.086*sab(k)/ab(k)
                else
                  sabo(k) = 2.0
                endif  
              endif
              if (ab(k).gt.0) then 
                 if (magtyp(1:1).eq.'A') abo(k)=-2.5*dlog10(ab(k))-48.60
                 if (magtyp(1:1).eq.'V') abo(k)=-2.5*dlog10(ab(k))+zp(k)
              else
                abo(k)=99.0
              endif
            else
              abo(k) = ab(k)
              sabo(k)= sab(k)
            endif
         enddo    
         if (cattyp(1:4).eq.'LONG') then
            j=2*imagm+2
            if (gbcont.eq.-1) then 
               read(paravc(j),'(i10)') cont
            else
               cont=gbcont
            endif   
            j=j+1
            call check_float(paravc(j),str_ch)
            read(str_ch,'(E20.9)') zs
            j=j+1
            str_inp = ' '
            do k = j,test
              str_inp = str_inp(1:lnblnk(str_inp)) // " " //
     >                 paravc(k)(1:lnblnk(paravc(k)))
             enddo
c            write(UO,'(I10,1x,A)') lnblnk(str_inp),str_inp 
         else
           if (gbcont.eq.-1) then 
              cont = 0
           else
              cont = gbcont
           endif   
         endif   
c         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  read obs mag catalogue 
      nobj = nobj + 1
      pass = 0
      abs_mag=-99
      redoing=0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  CONTEXT translated in used bands bused(k)=0 (NOT USED) or 1 (USED)
      nbused=0
      new_cont=0
      nbus=0
      nbul=0
      do k = 1, imagm
        if (cont.eq.0) then   ! by default if it is not read
           bused(k)=1
           nbused=nbused+1
           if (sab(k).le.0) nbul=nbul+1
        else
           bused(k)=bdincl(k-1,cont,imagm-1)
c  remove predefined context with forbitten bands 
           if(contforb.gt.0.and.
     .        bdincl(k-1,contforb,imagm-1).eq.1) bused(k)=0
           if (bused(k).eq.1) new_cont=new_cont+ 2**(k-1)
c
           if (bused(k).eq.1) nbused=nbused+1
           if (bused(k).eq.1 .and. sab(k).gt.0) nbus=nbus+1
           if (bused(k).eq.1 .and. sab(k).le.0) nbul=nbul+1
        endif
c same for Band used for scaling 
        if (bdscal.eq.0) then 
           buscal(k) = 1
        else
           buscal(k)=bdincl(k-1,bdscal,imagm-1)
        endif
      enddo   
c      Mag abs context 
      do k = 1,imagm
        macont=magabscont(k)
        do j = 1,imagm 
          mbused(k,j)=0
c          mbused(k,j)=bused(k)
c          if(macont.gt.0.and.
c     >     bdincl(j-1,macont,imagm-1).eq.0) mbused(k,j)=0
          if (macont.le.0)  mbused(k,j)=bused(j)
          if (macont.gt.0.and.
     >       bdincl(j-1,macont,imagm-1).eq.1) mbused(k,j)=1
        enddo
      enddo   
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Conversion from observed magnitude  to flux if it is the case  
c     Vega:  mi = -2.5 logfi + ZPi , AB:  mi = -2.5 logfi  - 48.60
c     fi = 10**-0.4(mi+zpi)     df =  f * dm/1.086
c     Re-scaling of errors in quadratic
c 
      if (typm.eq.'M' .or. typm.eq.'m') then 
         do k = 1,imagm
            if (magtyp(1:1).eq.'A') ab(k)=10**(-0.4*(ab(k)+48.60))
            if (magtyp(1:1).eq.'V') ab(k)=10**(-0.4*(ab(k)-zp(k)))
            if (sab(k).gt.0.0) then 
              if (min_err(k).gt.0) then 
                 sab(k) = dsqrt(sab(k)**2+min_err(k)**2)
              endif   
              sab(k) = ab(k)*sab(k)/1.086*fac_err
            endif
         enddo   
      else
         do k = 1,imagm
           if ( min_err(k).gt.0) then
             if (ab(k).gt.0.and.sab(k).gt.0.and.sab(k).ne.ab(k)) then
               sab(k)=1.086*sab(k)/ab(k)
               sab(k)=dsqrt(sab(k)**2+min_err(k)**2)
               sab(k)=ab(k)*sab(k)/1.086*fac_err
             endif   
           endif   
c           write(*,*) ab(k),sab(k)
         enddo      
      endif
c
ccccccccccAUTO-ADAPT
      if ((autoadapt(1:1).eq.'Y'.or.autoadapt(1:1).eq.'y').and.
     >     realise.lt.1) then
c        auto-adapt : keep only bright galaxies 
         if (ab(fl_auto).le.0.d0) goto 1
         if (-2.5*dlog10(ab(fl_auto))-48.60.ge.auto_thresmax) goto 1
         if (-2.5*dlog10(ab(fl_auto))-48.60.le.auto_thresmin) goto 1
c        auto-adapt : remove objets in upper-limit in filter adapt
         if (sab(fl_auto).le.0.d0.or.bused(fl_auto).eq.0) goto 1
c        If training also on color
         if (degre.gt.1) then
          if (sab(fl1).le.0.d0.or.bused(fl2).eq.0) goto 1
          if (sab(fl1).le.0.d0.or.bused(fl2).eq.0) goto 1
         endif
      endif
ccccccccccAUTO-ADAPT

ccccccccccTRY C17 METHOD
c     Compute the observed color in mag and the associated error
      if(nbC17.eq.imagm)then 
       do k = 1, imagm
        if (bused(k).eq.1  .and. sab(k).gt.0) then
         colC17(k) = -2.5*log10( ab(k)  / ab(C17fl(k)) ) 
         errC17(k) = sqrt((sab(k)/ab(k))**2.
     .                   +(sab(C17fl(k))/ab(C17fl(k)))**2.)
        endif
       enddo

      endif
ccccccccccTRY C17 METHOD
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c   Minimizes the search box in color space 
c   Extract record to be used according to colors 
      if (fastmod(1:1) .eq. "Y" .or. fastmod(1:1).eq."y") then
        call sear_rec(ab,sab,imagm,bused,sigcol,numcol,sel,
     >              fsort,f_index,colibmax,reclist,numrec)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccc
 36   if (redoing.eq.1)   liblength=colibmax
      mod_dist = 0
      do k = 1,100         ! Number of secondary peaks
        pdz(k) = 0
        nb(k)  = 1 
      enddo
      do k = 1, 50
         dzpdz(k) =0
      enddo
      zbay=-99.
      zbayi=-99.
      zbays=-99.
      do k = 1,imagm       ! Number of filters 
        goodfilter(k) = -1
        kap(k)    = -999
        mabs(k)   = -999
        kapq(k)   = -999
        mabsq(k)  = -999
        magm(k)   = -999
      enddo 
      if (outsp(1:1).eq.'Y'.or.outsp(1:1).eq.'y') then    ! Init. spectra
        do k = 1, wmax
          wsp(k)=0.
          fsp(k)=0.
          wq(k) =0.
          fq(k) =0.
          wst(k)=0.
          fst(k)=0.
          do j = 1,5 
            fgal(j,k)=0.
            wgal(j,k)=0.
          enddo
        enddo
      endif
      do k = 1,chimax    !   chi with z = 0 --> zmax 
         chi(2,k)= 1.e9
         maxlz(k)= 0
         recb(k) = -999
         imasb(k)= -999
         extb(k) = -99.
         ageb(k) = -99.
         zb(k)   = -99.
         dmb(k)  = -99.
         zfb(k)  = -99.
         mag_absb(k)=-99.
         chibay(k)=0. 
         do j = 1,imagm
           kcorb(j,k)=-9999.
         enddo
      enddo
c
      do nlib = 1 , 3     !  best fittting for the 3 libraries
        chimin(nlib)   = 1.e10
        zmin(nlib)     = -99
        dmmin(nlib)    = -99
        recmin(nlib)   = -99
        imasmin(nlib)  = -999
        agemin(nlib)   = -99.
        extmin(nlib)   = -99.
        zfmin(nlib)    = -99.
        mag_abs(nlib)  = -99.
        mod_distb(nlib)= -99.
      enddo
      if (fastmod(1:1) .eq. "Y" .or. fastmod(1:1).eq."y") then
         liblength=numrec
      else
         liblength=colibmax
      endif
      if (redoing.eq.1)   liblength=colibmax
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c START BIG LOOP 
      do i = 1, liblength    ! library models (correspondance with goto 2 )
         if ((fastmod(1:1) .eq. "Y" .or. fastmod(1:1).eq."y") 
     >       .AND. redoing.eq.0) then         
            model= modlib(reclist(i))
            exti = extilib(reclist(i)) 
            z    = zlib(reclist(i))
            age  = agelib(reclist(i))
            zfmod= zflib(reclist(i))
            do k = 1, imagm
              mag(k) = maglibf(k,reclist(i))
              magm(k) = maglib(k,reclist(i))
              kcor(k) = klib(k,reclist(i))
            enddo          
         else
           model = modlib(i)
           exti = extilib(i) 
           z    = zlib(i)
           age  = agelib(i)
           zfmod= zflib(i)
           do k = 1, imagm
             mag(k)  = maglibf(k,i)
             magm(k) = maglib(k,i)
             kcor(k) = klib(k,i)
           enddo   
         endif   

cccccccccccccccccccccccccccccccccccccccccccccc
c       Add if you don't want to compute photometric if zs wasn't mesured
c        if ( (zfix(1:1).eq.'Y'.or. zfix(1:1).eq.'y')
c     >       .and.(zs.le.0.001 .or. zs.ge.zmax))goto 1
cccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccAUTO-ADAPT
c       auto-adapt : keep object in a given model bin (ADAPT_MODBIN) 
c                            and in a given redshift bin (ADAPT_ZBIN)
        if ( (autoadapt(1:1).eq.'Y'.or. autoadapt(1:1).eq.'y')
     >    .and.(zs.lt.adzmin .or. zs.gt.adzmax .or. zs.ge.zmax)
     >    .and.realise.le.1) goto 1
ccccccccccAUTO-ADAPT

c 
        if (model.ge.0 .and. z.le.1.e-5) irecz0=i     ! for ABS Mag at z=0
c
c   do not read library  if  0<Zfix<9  and (Zlib-Zspec) > zstep/2 
        if ( (zfix(1:1).eq.'Y'.or. zfix(1:1).eq.'y')
     >    .and. zs.ge.0 .and. zs.lt.zmax 
     >    .and. model .lt. 1000 .and. model .gt.0) then
             if ( SNGL(dabs(z-zs)).gt.(SNGL(zstep/2.)) ) goto 2
        endif
c
c   distance modulus 
        do k = 1, chimax
            if (DABS(z-chi(1,k)).le.1.e-9) then
               mod_dist=dist_mod(k)
               goto 5
            endif
         enddo   
c  do not library beyond zmax 
 5       if (model.ge.0 .and. z.gt.zmax) goto 2
c
ccccccccccccccccccccc 
c  Defined the type of objects nlib=1 gal ; nlib=2 qso ; nlib=3 star
         if (model.lt.0)                    nlib = 3
         if (model.gt.1000)                 nlib = 2
         if (model.ge.1 .and. model.le.999) nlib = 1
c   do not use galaxy library  if  zgal > zmax_gal
         if (nlib.eq.1 .and. z.gt.zmax_gal ) goto 2 


ccccccccccAUTO-ADAPT
cccccccccccccccccccccccccccccccccccccccccc
c   auto-adapt : fix z , only gal librairy
        if ( (autoadapt(1:1).eq.'Y'.or. autoadapt(1:1).eq.'y')
     >       .and.realise.le.1) then
             if(sngl(dabs(z-zs)).gt.sngl((zstep/2.))) goto 2
             if(nlib.ne.1) goto 2
        endif
ccccccccccAUTO-ADAPT

c
ccccccccccccccccccccccccccccccccccccccccc
c   Measurement of scaling factor dm only with (fobs>flim), dchi2/ddm = 0
c   Uses only bands where observed fluxes is > 0 
         avmago = 0. 
         avmagt = 0. 
         nf = 0 
         do k = 1,imagm
c           if (ab(k).gt.0.and.sab(k).gt.0.and.
c     >         bused(k).eq.1 .and.buscal(k).eq.1) then
           if (ab(k).gt.0.and.sab(k).gt.0.and.
     >         bused(k).eq.1) then
             nf = nf + 1
             if(buscal(k).eq.1)then
               avmago = avmago + ab(k)*mag(k)/sab(k)**2
               avmagt = avmagt + mag(k)*mag(k)/sab(k)**2
             endif
           endif
         enddo
         if (nf.eq.0 .and. i.eq.1) then
             write(UO,*) spec,'  WARNING: No scaling --> No z' 
             write(UO,'(I10,1x,20(f6.2,1x))') spec,(sab(k),k=1,imagm)
             goto 4          ! go to write output 
         elseif (nf.eq.1 .and. i.eq.1 ) then 
             write(UO,*) spec,'  WARNING: 1 band only --> No z'
             goto 4 
         else   
            if (avmagt.lt.1.e-60) goto 2
             dm = avmago/avmagt
         endif
c
cccccccccccccccccccccccccccccccccccccccccc
c  Use upper-limits  for model rejection  
         do k = 1,imagm 
            if (bused(k).eq.1.and.sab(k).le.0.and.
     >           (dm*mag(k)).gt.ab(k)) goto 2
         enddo
cccccccccccccccccccccccccccccccccccccccccc
c  Use Mass scaling for model rejection
         if (nlib.eq.1 .and. lmasi.gt.0 .and. lmass.gt.0 ) then
           if (dm.lt.lmasi.or.dm.gt.lmass) goto 2
         endif
cccccccccccccccccccccccccccccccccccccccccc
c  Use Age at z larger than Age for Zform
         zform2 = zform(model)
         if (nlib.eq.1 .and. zform2.ge.z .and. age.gt.1.e4) then 
              tuniv  = timy(z,h0,om0,l0)      ! zinf -> z
              tzform = timy(zform2,h0,om0,l0) ! zinf -> zform
              tused  = tuniv-tzform           ! Age gal at z for a zform
              if (age.le.tused) goto 2
         endif   
ccccccccccccccccccccccccccccccccccccccccc
c   Abs Magnitude from model @ z=0 for rejection if babs defined 
         if (nlib.le.2 .and. babs.gt.0) then
c            abs_mag=magm(babs)-mod_dist-kcor(babs)
c            abs_mag= abs_mag-2.5*dlog10(dm)
            if (fastmod(1:1).eq."Y" .or. fastmod(1:1).eq."y" ) then 
               abs_mag = maglib(babs,reclist(irecz0)) -2.5*dlog10(dm)
             else
               abs_mag = maglib(babs,irecz0) -2.5*dlog10(dm)
             endif  
ccccccccccccccccccccccccccccccccccccccccc
c  rejection for galaxy based on Mag abs 
           if (nlib.le.1.and. 
     >         magabsl(1).lt.0 .and. magabsl(2).lt.0) then 
c     >   .and. magm(babs).le.100 .and. DABS(kcor(babs)).le.100) then
             if (magabsl(1).le.magabsl(2)) then
                if (abs_mag.lt.magabsl(1).OR.
     >              abs_mag.ge.magabsl(2)) goto 2
             elseif (magabsl(1).gt.magabsl(2)) then
                if (abs_mag.lt.magabsl(2).OR.
     >              abs_mag.ge.magabsl(1)) goto 2
             endif
           endif
         endif  
ccccccccccccccccccccccccccccccccccccccccc
c  Chi2 weighted by Nz if prior used
c         if (nlib.eq.1 .and. bp.gt.0 .and. bused(bp).eq.1 .and.
c     >       bused(bp_B).eq.1.and. bused(bp_I).eq.1 ) then
          if (nlib.eq.1 .and. bp.gt.0 .and. bused(bp).eq.1 
     >        .and.sab(bp).gt.0) then
           if (magtyp(1:1).eq.'A') then 
              iab = -2.5*dlog10(ab(bp)) - 48.60
c              color_rf=maglib(bp_B,irecz0)-maglib(bp_I,irecz0)
           elseif (magtyp(1:1).eq.'V') then 
               iab = -2.5*dlog10(ab(bp))+zp(bp) + abcor(bp) 
c               color_rf=maglib(bp_B,irecz0)-maglib(bp_I,irecz0) 
c     >                 +abcor(bp_B)-abcor(bp_I)
           endif
             pweight=nzpriorVVDS5(iab,model,z)
c             pweight=nzprior2(iab,color_rf,z)
c           if ( gallib(1:lnblnk(gallib)) .eq. 'LIB_AVEROI' ) then
c             pweight = nzprior4(iab,model,z)
c           else
c             pweight=nzprior3(iab,z)
c           endif
         endif

         if(nbC17.ne.imagm)then !TRY C17 METHOD
ccccccccccccccccccccccccccccccccccccccccc
c Measurement of chi^2
           chi2 = 0.
           do k = 1, imagm
             if (bused(k).eq.1  .and. sab(k).gt.0) then
               chi2 = chi2 + ((ab(k)-dm*mag(k))/sab(k))**2
             endif
           enddo
c
ccccccccccTRY C17 METHOD 
c        Compute chi2 with color according to the color terms defined
c        in C17_METHOD
         else
           chi2 = 0.d0
           do k = 1, imagm
            if(bused(k).eq.1.and.sab(k).gt.0.and.
     >      errC17(k).gt.0.and.C17fl(k).gt.0)then
              chi2=chi2 + 
     >  ((colC17(k) + 2.5*log10(mag(k) / mag(C17fl(k)))) / errC17(k))**2
            endif
           enddo
         endif
ccccccccccTRY C17 METHOD
c
         if (nlib.eq.1 .and. bp.gt.0 ) then            
             if(pweight.le.0.d0)then 
                chi2 = 1.d10
             else         
                chi2 = chi2 - 2*dlog(pweight)
                if(chi2.lt.0.d0)then
c                  write(*,*)"ALARM CHI2"
c                  write(*,*)chi2,spec,z
                  chi2=0.d0
                endif
             endif
         endif     
cccccccccccccccccccccccccccccccccccccccccc
c Extimate of the best chi at any redshift for the 3 types of objects
         if (chi2.lt.chimin(nlib) ) then 
            chimin(nlib) = chi2
            zmin(nlib)   = z
            dmmin(nlib)  = dm
            recmin(nlib) = i
            recmin0(nlib)= irecz0  ! to keep magnitude at z=0
            imasmin(nlib)= model
            agemin(nlib) = age
            extmin(nlib) = exti
            zfmin(nlib)  = zfmod
            mag_abs(nlib)=abs_mag
            if (nlib.le.2) then 
               mod_distb(nlib) = mod_dist
            endif   
         endif
 3       continue
ccccccccccccccccccccccccccccccccccccccccc
c study of the global chi2 evolution for gal library only 
         if (chi2.lt.1.e9.and.nlib.eq.1) then
           do nchi = 1,chimax
              if (DABS(z-chi(1,nchi)).le.1.e-5) then
c baysian 
c      write(*,*) 'OK',nobj,z,nchi,chi2 
c      write(*,*) '  ' 
                  chibay(nchi) = chibay(nchi)+dexp(-0.5*chi2)
c            write(*,*) '  ',nobj,i,z,chi2,chibay(nchi)
                  if (chi2.lt.chi(2,nchi)) then 
                     chi(2,nchi) = chi2
                     recb(nchi)  = i
                     recb0(nchi) = irecz0  ! to keep magnitude at z=0
                     imasb(nchi) = model
                     extb(nchi)  = exti
                     do k = 1,imagm
                        kcorb(k,nchi)=kcor(k) 
                     enddo
                     zfb(nchi)   = zfmod
                     ageb(nchi)  = age
                     zb(nchi)    = z
                     dmb(nchi)   = dm
                     mag_absb(nchi)= abs_mag 
                  endif
                  goto 2
              endif
           enddo
         endif
 2       continue    
      enddo                 ! END of LOOP for Chi^2 measurement 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Analysis of ML function .  Normalizing to peak at ML=1
      ndz=1
      if (chimin(1).ge.1e9 .or.
     >  ((zfix(1:1).eq.'Y' .or. zfix(1:1).eq.'y')
     >  .and. zs.ge.0 .and. zs.lt.zmax .and. model.lt.1000) ) then
        goto 4 
      endif
      if (chimin(1).lt.1e9) then
        do k = 1,chimax
          xp(k) = chi(1,k)
          yp(k) = dexp(-0.5*(chi(2,k) - chimin(1)))
          if (yp(k).le.1.e-9) yp(k) = 0.d0
        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Redshift  uncertainties (zmin,zmax) for dChi2=1,2.71,6.63
        chibest=chimin(1)
        dchi = 1.0
        call DCHI2(chi,chimax,chibest,dchi,z68i,z68s)
        dchi = 2.71
        call DCHI2(chi,chimax,chibest,dchi,z90i,z90s)
        dchi = 6.63
c        dchi = 9.0
        call DCHI2(chi,chimax,chibest,dchi,z99i,z99s)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c   Extimates of the secondary peak in z
        call ZPEAK_SCALE(chi,yp,chimax,zstep,dz_win,min_thres,nb,ndz)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   normalizing Pdz curve to 1
        z0=0.
        call TRAPZD(xp,yp,chimax,z0,zmax,mlarea) 
c   normalizing baysian zphot curve to 1
        call TRAPZD(xp,chibay,chimax,z0,zmax,barea)
c        write(*,*) '   '
c        write(*,*) nobj,' area=',barea
c        write(*,*) '   '
        if (barea.gt.0) then 
c   redshift uncertainties at 68%  for baysian zphot
          call PROBAZBAY(xp,chibay,chimax,barea,zmed,zpdzi,zpdzs)
          zbay = zmed
          zbayi=zpdzi(2)     
          zbays=zpdzs(2)     
        endif
c        write(*,*) zbay,zbayi,zbays
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Redshift uncertainties (zmin,zmax) for 68%,90%,99%
        call PROBAZ(xp,yp,chimax,mlarea,zpdzi,zpdzs,probz)
        zml68i=zpdzi(1)
        zml68s=zpdzs(1)
        zml90i=zpdzi(2)
        zml90s=zpdzs(2)
        zml99i=zpdzi(3)
        zml99s=zpdzs(3)
c
        do k = 1, chimax
          yp(k) = yp(k) / mlarea
        enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Computing fraction included between +/- Dz = 0.1*(1+z)
        dzp=0.1
        dzml = dzp*(1+zmin(1))
        zinf = zmin(1) - dzml
        zsup = zmin(1) + dzml
        if (zinf.lt.0)    zinf = 0
        if (zsup.gt.zmax) zsup = zmax
        call TRAPZD(xp,yp,chimax,zinf,zsup,summl)
        pdz(1) = summl*100.  
        if (ndz.gt.1) then 
          do k = 2 , ndz
             dzml = dzp*(1+zb(nb(k)))
             zinf = zb(nb(k)) - dzml
             zsup = zb(nb(k)) + dzml
             if (zinf.lt.0)    zinf = 0
             if (zsup.gt.zmax) zsup = zmax
             call TRAPZD(xp,yp,chimax,zinf,zsup,summl)
             pdz(k) = summl*100.                  
          enddo
        endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Computes proba in various intervals of z
        if (npdz.gt.1) then
           j = 1
           do k = 1, npdz, 2
              zinf=int_pdz(k)
              zsup=int_pdz(k+1)
              call TRAPZD(xp,yp,chimax,zinf,zsup,summl)
              dzpdz(j)=summl*100
              j = j + 1
           enddo
        endif    
        do k = 1, chimax
          yp(k) = yp(k) * mlarea
        enddo
      else
         pdz(1) = -99
        if (npdz.gt.1) then
           j = 1
           do k = 1, npdz, 2
              if (int_pdz(k).le.zmin(1).and.
     >             int_pdz(k+1).gt.zmin(1)) then
                 dzpdz(j) = 100
              else
                  dzpdz(j) = 0
              endif
              j=j+1
           enddo
        endif    
      endif
c

c
c Test if need to be redone 
 4    if ( (fastmod(1:1).eq."Y" .or. fastmod(1:1).eq."y") 
     >   .AND. chimin(1).gt.1e8  .AND. redoing.eq.0 ) then
        redoing=1
        goto 36
      endif  
ccccccccccccccccccccccccccccccccccccccc
c  writing output   (correspondance with goto 4)
c
c  Converting fluxes to Magnitude   
      do k = 1, imagm
         if ( ab(k).gt.0 .and. sab(k).gt.0 ) then 
           sab(k) = 1.086*sab(k)/ab(k)
           if (magtyp(1:1).eq.'A') then 
              ab(k)=-2.5*dlog10(ab(k))-48.60
           elseif (magtyp(1:1).eq.'V') then 
              ab(k)=-2.5*dlog10(ab(k))+zp(k)
           endif   
         elseif (ab(k).gt.0 .and. sab(k).le.0 ) then
             sab(k) = -1
             if (magtyp(1:1).eq.'A') then
                ab(k)=-2.5*dlog10(ab(k))-48.60
             elseif (magtyp(1:1).eq.'V') then
                ab(k)=-2.5*dlog10(ab(k))+zp(k)
             endif   
         elseif (ab(k).le.0 ) then
            sab(k) = 99.
            ab(k) = 99.
         endif
      enddo
cccccccccccccccccccccccccccccccccccccccccc
c  Estimation of the kappa correction kap(mod) = mabs(z) - mabs(0)
c        kap = mapp(z) - fz(z) - (mapp(0) -fz(0))  
c   and apparent magnitudes in z and z = 0 
      if (zmin(1).ge.0) then
         do k = 1, imagm
            if (fastmod(1:1).eq."Y" .or.fastmod(1:1).eq."y" ) then
               magm0(k) = maglib(k,reclist(recmin0(1)))
             else    
               magm0(k) = maglib(k,recmin0(1))
             endif 
         enddo     
c  Interpolation of Z_best (zmin(1)) via Chi2 curves
         if (zintp(1:1).eq."y".or.zintp(1:1).eq."Y") then
             zbest=zmin(1)
             call int_parab(chi,chisize,chimax,zbest,zintb)
             zmin(1)=zintb
             mod_distb(1) = funz(zintb,h0,om0,l0)
         endif
c   Z_best replaced by Z_spec if defined 
         if( (zfix(1:1).eq.'Y' .or. zfix(1:1).eq.'y')
     >.and. zs.ge.0 .and. zs.lt.zmax .and. model.lt.1000) then 
             zmin(1)=zs
             mod_distb(1) = funz(zs,h0,om0,l0)
         endif
c   Interpolation for k-correction  and theoretical magnitudes if needed
         if(fastmod(1:1).eq."Y" .and.fastmod(1:1).eq."y" )then
             index=reclist(recmin(1))
         else
             index=recmin(1)
         endif
            call kInterp(zmin(1),index,imagm,zlib,modlib,klib,kap)
            call kInterp(zmin(1),index,imagm,zlib,modlib,maglib,magm)
c  Computes  ABSOLUTE MAGNITUDES  FOR GALAXIES with various METHOD 
c         write(*,*) (method(k),k=1,imagm)
         call absMagPro(method,spec,imagm,bapp,zstep,dz,h0,om0,l0,
     >            fobs,fobs4,mbused,zmin(1),dmmin(1),ab,magm,sab,kap,
     >            magm0,minkcol,goodfilter,minkcolor,mabs)
c         write(*,*) (method(k),k=1,imagm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
         do k = 1, imagm 
            magm(k) = magm(k) -2.5*dlog10(dmmin(1))
         enddo   
         
      else                 ! if (Chi^2>1.e9 & no Zphot) 
         do k = 1, imagm
            kap(k)  = -99.99
            mabs(k) = -99.99
            magm(k) = -99.99
         enddo   
      endif
cccccccccccccccccccccccccccccccccccccccccc
c Mag_abs for QSO 
      if (zmin(2).ge.0) then 
         do k = 1, imagm
            if (fastmod(1:1).eq."Y" .or.fastmod(1:1).eq."y" ) then
               kapq(k)   = klib(k,reclist(recmin(2)))
c               magm0(k)  = maglib(k,reclist(recmin0(2)))
             else    
               kapq(k)   = klib(k,recmin(2))
c               magm0(k)  = maglib(k,recmin0(2))
             endif    
c  Compute :  MAG_ABS derived from MAG_OBS in same filter. 
            mabsq(k)=ab(k)-mod_distb(2)-kapq(k)
         enddo              
      else
        do k = 1, imagm
            kapq(k)  = -99.99
            mabsq(k) = -99.99
         enddo   
      endif

ccccccccccAUTO-ADAPT
c     auto-adapt : keep information on pedicted apparent magnitudes  
      if (zmin(1).ge.0.d0) then
       if ((autoadapt(1:1).eq.'Y'.or. autoadapt(1:1).eq.'y')
     .      .and.realise.le.1)then
        if(ngals_ada.lt.zadapt)then
         cont_ada(ngals_ada)=new_cont
         if(imasmin(1).lt.admmin .or.imasmin(1).gt.admmax)goto 1
         mod_ada(ngals_ada)=imasmin(1)
         zs_ada(ngals_ada)=zs
         do k=1,imagm
          ab_ada(ngals_ada,k)=ab(k)
          sab_ada(ngals_ada,k)=sab(k)
          magm_ada(ngals_ada,k)=magm(k)
         enddo
        endif
        ngals_ada=ngals_ada+1
        endif
      endif
ccccccccccAUTO-ADAPT

10001 format(i9,1x,f9.3,1x,i6,1x,E12.6,1x,f9.3,1x,70(f12.5,1x))

c
cccccccccccccccccccccccccccccccccccccccccccc
c objects measured 
      do k = 1,3 
         if (chimin(k).lt.1e9) nmeas(k)=nmeas(k) + 1
      enddo  
cccccccccccccccccccccccccccccccccccccccccc
c writing in output file 
      if(realise.ge.2.or.(autoadapt(1:1).eq.'N'.or.
     >   autoadapt(1:1).eq.'n'))then
        call WRITE_OUT(wpara,iwpara,spec,
     >  zmin,chimin,imasmin,agemin,extmin,zfmin,dmmin,pdz,mag_abs,
     >  cont,zs,str_inp,zb,
     >  chi,imasb,ageb,extb,zfb,dmb,mag_absb,nb,
     >  abo,sabo,kap,mabs,imagm,goodfilter,magm,
     >  dzpdz,npdz,
     >  z68i,z68s,z90i,z90s,z99i,z99s,
     >  zml68i,zml68s,zml90i,zml90s,zml99i,zml99s, 
     >  zbay,zbayi,zbays,
     >  mabsq,kapq, 
     >  nbus,nbul,
     >  mod_distb,
     >  str_out,iwout)
c
        write(4,'(200(A,2x))') 
     >  (str_out(k)(1:lnblnk(str_out(k))),k=1,iwout)
ccccccccccccccccccccccccccccccccccccccccc
c writing PDZ files 
        if (outpdz(1:4).ne."NONE" .AND. outpdz(1:4).ne."none"
     >     .and. fastmod(1:1).ne."Y" .AND.fastmod(1:1).ne."y") then
          do k = 1,chimax
             if (yp(k).le.1.e-10) then
                do i = 1, nfabs
                   pdz_mabs(i,k) =-99.99
                enddo
             else 
                pdz_z = chi(1,k)
                pdz_dm= dmb(k)
                pdz_distm=dist_mod(k)
                do i = 1,imagm 
                   pdz_kcor(i)= kcorb(i,k)
                   magm0(i)   = maglib(i,recb0(k))
                   magm(i)    = maglib(i,recb(k))
                enddo
c  Computes  ABSOLUTE MAGNITUDES  FOR GALAXIES with various METHOD 
                do i = 1,nfabs
                 fabs=pdz_fabs(i)
                 call absMag1F(method,spec,fabs,bapp,zstep,dz,
     >      pdz_distm,fobs,fobs4,mbused,pdz_z,pdz_dm,ab,magm,sab,
     >      pdz_kcor,magm0,minkcol,goodfilter,minkcolor,pdz_mabsz)
                 pdz_mabs(i,k) = pdz_mabsz
               enddo
             endif
          enddo
          write(44,'(500(f8.3,1x))') (yp(k),k=1,chimax)
          write(45,'(500(i6,1x))')   (imasb(k),k=1,chimax)
          do i = 1,nfabs
            ufabs=45+i 
            write(ufabs,'(500(f8.3,1x))') (pdz_mabs(i,k),k=1,chimax)
          enddo
        endif
      endif

cccccccccccccccccccccccccccccccccccccccccc
c info in screen
      if (UO .eq. 6) then    
      if(autoadapt(1:1).eq.'N' .or. autoadapt(1:1).eq.'n'
     >    .or. realise.ge.2 )then
      if (fastmod(1:1).eq."Y" .or. fastmod(1:1).eq."y") then
         write(UO,603) spec,numrec,zmin(1),pdz(1),mag_abs(1),
     >   zmin(2),mag_abs(2),imasmin(3),char(13)
         call flush(UO)
      else
         write(UO,602) spec,zmin(1),pdz(1),mag_abs(1),
     >   zmin(2),mag_abs(2),imasmin(3),char(13)
         call flush(UO)
      endif   
      endif
      endif
c
 602  format("object ->",I10,1x,5(f6.2,1x),I4,a1,$)
 603  format("object ->",2(I10,1x),5(f6.2,1x),I4,a1,$)
cccccccccccccccccccccccccccccccccccccc
c  Writing the spectra and chi2 curves  
      if (outsp(1:1).eq.'Y'.or.outsp(1:1).eq.'y') then
c output file
        write(ospec,'(i9.9)') spec 
        ospec= 'Id' // ospec(1:lnblnk(ospec)) // '.spec'
        lbdminwr = 1.e20
        lbdmaxwr = 0 
        do i = 1, imagm 
          if ( (flmoy(i)-flwidth(i)) .le. lbdminwr ) then
            lbdminwr = flmoy(i) - flwidth(i)
          endif
          if ( (flmoy(i)+flwidth(i)) .ge. lbdmaxwr ) then
            lbdmaxwr = flmoy(i) + flwidth(i)
          endif
        enddo
c GALAXY 
        nspmaxg=0
        nspmaxq=0
        nspmaxs=0
        nspmax =0
        if (chimin(1).lt.1e9 .and. zmin(1).ge.0) then
c           write(*,*) 'read galaxy spectrum ...'
c
          call read_spec(gallib,zmin(1),imasmin(1),agemin(1),
     >    extic,iext,extmin(1),opal,opat,iopa,
     >    lbdminwr,lbdmaxwr,1,wsp,fsp,nsp)
c
          ztemp=zmin(1)
          mod_dist2=funz(ztemp,h0,om0,l0)
c
          do i = 1, nsp
             fgal(1,i)= fsp(i)-2.5*dlog10(dmmin(1))
     >                       + mod_dist2
             wgal(1,i)= wsp(i)
          enddo
          if (nsp.gt.nspmaxg )   nspmaxg=nsp
          if (ndz.gt.1) then
             if (ndz.gt.5) ndz=5   ! just take the 5 first spectra
            do k = 2, ndz  
              call read_spec(gallib,zb(nb(k)),imasb(nb(k)),ageb(nb(k)),
     >          extic,iext,extb(nb(k)),opal,opat,iopa,
     >          lbdminwr,lbdmaxwr,2,wsp,fsp,nsp) 
              ztemp=zb(nb(k))
              mod_dist2=funz(ztemp,h0,om0,l0)
              do i = 1, nsp
                fgal(k,i)= fsp(i)-2.5*dlog10(dmb(nb(k)))
     >                       +mod_dist2
                wgal(k,i)= wsp(i)
              enddo  
              if (nsp.gt.nspmaxg)    nspmaxg=nsp
            enddo
          endif
          if (nspmaxg.gt.nspmax) nspmax=nspmaxg
        endif
c
c QSO
        if (chimin(2).lt.1e9 .and. zmin(2).ge.0) then
           call read_spec(qsolib,zmin(2),imasmin(2),agemin(2),
     >      extic,iext,extmin(2),opal,opat,iopa,
     >      lbdminwr,lbdmaxwr,4,wsp,fsp,nsp)
            ztemp=zmin(2)
            mod_dist2=funz(ztemp,h0,om0,l0)
            do i = 1, nsp
              fq(i)= fsp(i)-2.5*dlog10(dmmin(2))+mod_dist2
              wq(i)= wsp(i)
            enddo
          if (nsp.ge.nspmaxq) nspmaxq=nsp
          if (nspmaxq.gt.nspmax) nspmax=nspmaxq
        endif
c
c STAR
        if (chimin(3).lt.1e9 ) then
           call read_star(starlib,imasmin(3),lbdminwr,
     >                    lbdmaxwr,5,wsp,fsp,nsp)
            do i = 1, nsp
              fst(i)= fsp(i)-2.5*dlog10(dmmin(3)) 
              wst(i)= wsp(i)
            enddo
          if (nsp.ge.nspmaxs ) nspmaxs=nsp
          if (nspmaxs.gt.nspmax) nspmax=nspmaxs
        endif
c OUTPUT SPECTRA
        open(7,file=ospec,status='unknown')
        write(7,'(i6,1x,2(f8.3,1x),E12.6)') imagm,zs,zmin(1),chimin(1)
        do k = 1, imagm 
           if (magtyp(1:1).eq.'V') ab(k) = ab(k) + abcor(k)
          write(7,701) ab(k),sab(k),flmoy(k),flwidth(k),abcor(k)
        enddo
        if (ndz .eq. 1) then
          write(7,702) ndz,zmin(1),z68i,z68s,chimin(1),pdz(1),
     >        mag_abs(1),imasmin(1),agemin(1),extmin(1)
        else                      
         write(7,702) ndz,zmin(1),z68i,z68s,chimin(1),pdz(1),
     >         mag_abs(1),imasmin(1),agemin(1),extmin(1),
     >         (zb(nb(k)),chi(2,nb(k)),pdz(k),mag_absb(nb(k)),
     >         imasb(nb(k)),ageb(nb(k)),extb(nb(k)),k=2,ndz)
        endif
        write(7,703) zmin(2),chimin(2),mag_abs(2),imasmin(2),
     >                        imasmin(3),chimin(3)
c
        if (nspmax.gt.wmax) nspmax=wmax   ! must be lower than size declaration
        do i = 1, max(nspmax,chimax)
           if (i.gt.chimax) then
              val1 = 0.
              val2 = 0.
           else
              val1 = chi(1,i)
              if (chimin(1).lt.1.e9)   val2 = yp(i)
              if (chimin(1).ge.1.e9)   val2 = chi(2,i)
           endif   
           if (i.gt.nspmaxg) then
              do k = 1, ndz
                 wgal(k,i) = 0.
                 fgal(k,i) = 0.
              enddo  
           endif
           if (i.gt.nspmaxq) then
               wq(i) =0.
               fq(i) =0.
           endif
           if (i.gt.nspmaxs) then
               wst(i) =0.
               fst(i) =0.
           endif
           write(7,704) val1,val2,wq(i),fq(i),wst(i),fst(i),
     >                 (wgal(k,i),fgal(k,i),k=1,ndz) 
        enddo
        close(7)
      endif
c
 701  format(50(f10.3,1x))
 702  format(i6,1x,3(f9.3,1x),E12.6,2(1x,f9.3),1x,i8,1x,E12.6,1x,f9.3,
     > 1x,70(f9.3,1x,E12.6,2(1x,f9.3),1x,i8,1x,E12.6,1x,f9.3,1x))
 703  format(f9.3,1x,E12.6,1x,f9.3,2(1x,i6),1x,E12.6)
 704  format(f10.3,1x,E12.6,1x,16(E14.7,1x))

c
      


c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  LOOP ON THE NUMBER OF OBJECTS IN THE OBSERVED CATALOG 
 1    continue   ! iteration on the observ. catalogue 
      enddo
 8    nobjm = k - 1
      close(90)


ccccccccccAUTO-ADAPT
c
c     If finished (convergence in auto-adapt is reached or no auto-adapt) 
      if(autoadapt(1:1).eq.'N' .or. autoadapt(1:1).eq.'n'
     >    .or. realise.ge.2 )then

        write(UO,*)  '            '      
        write(UO,*)  ' Number of measured objects :'
        write(UO,'(4(A,I8,1x))') " GAL:",nmeas(1)," QSO:",nmeas(2),
     >"STAR:",nmeas(3)," over ",nobj 
        write(UO,*)  ' Results in file : ',outf(1:lnblnk(outf))
        write(UO,*)  ' That s all Folks !! '

      else

       ngals_ada=ngals_ada-1

c      keep free parameters before new step  
       do k=1,imagm
        a0in(k)=a0(k)
        a1in(k)=a1(k)
        a2in(k)=a2(k)
        a3in(k)=a3(k)
c        write(*,*) ab_ada(ngals_ada,k),sab_ada(ngals_ada,k)
       enddo

c      main procedure for auto-adapt
       write(41,*)"ITERATION",iteration
c       call auto_adapt(bused,chiin,chifit,a0in,a0,a1in,a1,
c     >             a2in,a2,a3in,a3,realise,degre,adcont,
c     >             residu,adapterror)
       call auto_adapt(chiin,chifit,a0in,a0,a1in,a1,
     >             a2in,a2,a3in,a3,realise,degre,adcont,
     >             residu,adapterror)
c      new chi2
       chiin=chifit

c      write at each step the corrections applied and the min error
c       if (realise.le.1) then
         do k = 1, imagm
           write(42,'(3(I4,1x),100(f12.5,1x))') realise,iteration,
     >     k,a0(k),a1(k),a2(k),a3(k),min_err(k)
         enddo
c       endif

c      keep the value for the smallest chi2
       if(residu.le.res_best)then
        res_best=residu
        iter_best=iteration
        do k=1,imagm
         a0best(k)=a0(k)
         a1best(k)=a1(k)
         a2best(k)=a2(k)
         a3best(k)=a3(k)
         min_errbest(k)=min_err(k)
        enddo
       endif
       

c      stop iteration
       if(realise.ge.2)then

c        write(UO,'(A)') "   " 
        write(UO,'(A)') "#####################################"
        write(UO,'(A)') " --> End of training : Corrections are:"
        write(UO,'(A20,1x,100(A9,1x))')
     >  "Filter","a0","a1","a2","a3"
           do k = 1, imagm
              write(UO,'(A20,1x,100(f9.3,1x))') 
     >          valf(k)(1:lnblnk(valf(k))),a0(k),a1(k),a2(k),a3(k) 
           enddo
c   writing coef in output 
           write(4,'(A)') "#####################################"
           write(4,'(A)') "# AUTO-ADAPT Coefficients:"
           write(4,'(A,1x,A20,1x,100(A9,1x))')
     >  "#"," Filter","a0","a1","a2","a3"
           do k = 1, imagm
              write(4,'(A,1x,A20,1x,100(f9.3,1x))') 
     >          "#",valf(k)(1:lnblnk(valf(k))),a0(k),a1(k),a2(k),a3(k) 
           enddo
           write(4,'(A,1x,100(f7.3,A))') 
     >            "# SHIFTS ",(a0(k),",",k=1,imagm)
           write(4,'(A)') "#####################################"
c
           write(UO,'(A)') "#####################################"
           write(UO,'(A)') " --> Running lephare ..."
       endif

       goto 6   !Redo the loop in auto-adapt

      endif
c
      if(autoadapt(1:1).eq.'Y'.or. autoadapt(1:1).eq.'y') then
         close(41)
         close(42)
      endif
      close(4)

      if (outpdz(1:4).ne."NONE" .and. outpdz(1:4).ne."none") then
        close(43)
        close(44)
        close(45)
        do i = 1,nfabs
          ufabs=45+i 
          close(ufabs)
        enddo
      endif

ccccccccccAUTO-ADAPT




      if (UO.ne.6)   close(UO)
     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  FORMATS 
 530  format(1x,f6.3,1x,E12.6,1x,f9.6,1x,E12.6)
 533  format(1x,I8,1x,i3,1x,4(f8.3,1x),2(1x,f11.6))
 534  format(1x,I4,10(1x,E13.6))
 535  format(3(1x,f6.3),1x,f10.3,1x,i8,1x,2(i3,1x),f9.5,1x,e12.6)
c 
c
      stop
      end
