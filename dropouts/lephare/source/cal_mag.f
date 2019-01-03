c     last modif 21/01/00
c     AUTHOR : S. ARNOUTS 
c
      subroutine cal_mag(flux,lbflux,nlf,rep,lbrep,nlr,
     >        flveg,lbveg,nlveg,mtyp,veg,ab,abveg)
c
c     This subroutine computes the mag in AB (ab) or Vega (veg) system
c      and gives the AB correction abveg as Mab= Mvega + abveg
c                for a given SED (flux,lbflux)1->nlf
c                and a given filter(rep, lbrep)1->nlr
      implicit none
      integer*4 wmax,maxsize
c      parameter (maxsize=110000,wmax=8000)
      parameter (maxsize=110000)
      INCLUDE 'dim_wave.decl'      
      integer*4 nlf,nlr,nlveg,j,p,smax,pmin
      real*8    flux(wmax),lbflux(wmax)
      real*8    lbveg(wmax),flveg(wmax)
      real*8    lbrep(wmax),rep(wmax)
c      real*8    ls(maxsize),fs(maxsize),ts(maxsize)
      real*8    ls(wmax),fs(wmax),ts(wmax)
      real*8    veg,ab,abveg
      real*8    fmel,arean,lmean,trans,fluxmean,area
      real*8    c,conv,dlbd,int_vega
      character mtyp*4096
      parameter(c=2.99792458e+18)      
c
       int_vega=0.
       pmin=0
       veg=0.d0
       ab=0.d0
       abveg=0.d0
       if (mtyp(1:1).eq.'V') pmin=1
       if (mtyp(1:1).eq.'A') pmin=2
       do p = pmin, 2
c compute Vega magnitude 
          fmel  = 0.d0
          arean = 0.d0
          area  = 0.d0     
c          write(*,*) nlveg,nlr 
          if (p.eq.1) then
         call sampling(lbveg,flveg,nlveg,lbrep,rep,nlr,ls,fs,ts,smax)  
          else
         call sampling(lbflux,flux,nlf,lbrep,rep,nlr,ls,fs,ts,smax)
          endif
c          smax = nlveg+nlr
c          do j = 1,smax
c            ls(j) = DBLE(j)+1000.d0
c            fs(j) = 1.d0
c            ts(j) = 1.d0
c          enddo
          do j = 1,smax - 1
             lmean = (ls(j) + ls(j+1))/2.
             trans = (ts(j) + ts(j+1))/2.
             fluxmean = (fs(j)+fs(j+1))/2.
             dlbd  =  ls(j+1)-ls(j) 
             conv= c/lmean**2        
             arean= arean + trans*dlbd*conv
             area = area  + trans*dlbd
             fmel  = fmel + fluxmean*trans*dlbd  
          enddo
c          write(*,*) area,arean,fmel
          if ( p .eq. 1) then
            if (fmel.gt.0.d0) then
              abveg = -2.5*dlog10(fmel/arean)-48.59
              int_vega = fmel
              veg = -2.5*dlog10(fmel/area)
            else
              abveg=999.99
              int_vega=0.d0
            endif
          elseif (p.eq.2 .and. pmin.eq.1) then
            if (fmel.gt.0) then
              ab  = -2.5*dlog10(fmel/arean) -48.59
              veg = -2.5*dlog10(fmel/int_vega)+0.03  
            else
              ab=999.99
              veg=999.99
            endif
          elseif (p.eq.2 .and. pmin.eq.2) then
             if (fmel.gt.0.d0) then
               ab = -2.5*dlog10(fmel/arean)-48.59
             else
               ab= 999.99
             endif
             veg=999.99 
             abveg=999.99
          endif      
       enddo
c
       return
       end
