c   last modif : 19/01/00
c     program read_star
      SUBROUTINE read_star(lib,modb,wavemin,wavemax,sampl,wave,fl,spmax)
c
c     Construction des differents modeles d'etoiles pour la 
c     mesure des couleurs d'etoiles et pour etre utiliser dans 
c     la mesure des z-photometriques 
c      Librairie stellaire vient de Pickles, 1998, PASP 110 
      implicit none
      integer*4 modb,nrec,recmin,recmax,recmax0,wmax
c      parameter (wmax=8000)
      INCLUDE 'dim_wave.decl'      
      integer*4 imas,i,j,k,iw,reclmax,nmod,nr
      integer*4 lnblnk,test,sampl,spmax,model
      real*8    wave(wmax),dummy
      real*8    fl(wmax)
      real*8    vec(wmax),vecl(wmax)
      real*8    wavemin,wavemax,c,val
      parameter(c=2.99792458e+18) 
      character*512  lib,file,zpdir,zpwork,paravc(100),param,str
c
      recmin=0
      recmax=0
      do k = 1, wmax 
        wave(k)= 0.
        fl(k)  = 0.
      enddo
      call getenv('LEPHAREDIR',zpdir)
      test=lnblnk(zpdir)
      if (test .eq. 0) then
        write(6,*) 'WARNING :  variable LEPHAREDIR not defined'
        stop
      endif
      call getenv('LEPHAREWORK',zpwork)
      test=lnblnk(zpwork)
      if (test .eq. 0) then
        write(6,*) 'WARNING :  variable LEPHAREWORK not defined'
        stop
      endif
c
      file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >       // lib(1:lnblnk(lib)) //'.doc'
c
      param='NUMBER_SED'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') imas
      param='RECORD_LENGTH'
      call read_para(param,file,paravc,test)
      read(paravc(1),'(i8)') reclmax      
      recmax0 = 0 
      do i = 1,imas
                  write(str,'(i6.6)') i 
         param = 'MOD_' // str(1:lnblnk(str))
         call read_para2(param,file,paravc,test)
         read(paravc(1),'(i8)') model
         if (test.ne.0) then
            read(paravc(3),'(i8)') nrec
           if (modb.eq.model) then
             recmin=recmax0+1
             recmax=recmin+ nrec
             goto 30
           else 
             recmax0 = recmax0 + nrec + 1
           endif
         else
            call err_option(param,1)
         endif   
      enddo   
c
 30   file = zpwork(1:lnblnk(zpwork)) //'/lib_bin/'
     >       // lib(1:lnblnk(lib)) //'.bin'
      open(27,file=file,status='unknown',access='direct',recl=reclmax)
      do j = recmin,recmax
        read(27,rec=j) nr,nmod,dummy,iw,(vec(k),k=1,iw)
        spmax = 0 
        if (dummy.eq.-1) then
          do k = 1, iw, sampl
             vecl(k)=vec(k)
             if (vec(k).gt.wavemin.and.vec(k).lt.wavemax) then
               spmax=spmax+1
               wave(spmax) = vec(k)
             endif
          enddo
        else
           val=99.0
          do k = 1, iw, sampl
             if (vecl(k).gt.wavemin.and.vecl(k).lt.wavemax) then
               spmax = spmax+1
               if (vec(k).le.0)  fl(spmax) = val
               if (vec(k).gt.0)  then
                  fl(spmax) = vec(k)*wave(spmax)*wave(spmax)/c
                  fl(spmax) = -2.5*dlog10(fl(spmax)) -48.59
               endif
             endif
          enddo             
        endif   
      enddo  
      close(27)
c
      return
      end




