cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine xweight(ngals_ada,ab_ada,sab_ada,weight)

      implicit none

      integer*4 i,nbf,total
      INCLUDE 'dim_filt.decl'
      integer*4 ngals_ada,nstep
      real*8 x,xcmp(100),xmin,xmax,xstep,weight(30000)
      real*8 ab_ada(30000,nbf),sab_ada(30000,nbf)


c     Calcule le weight a appliquer a chaque gal
      xmin=-10.d0
      xmax=10.d0
      xstep=0.25d0
      nstep=nint((xmax-xmin)/xstep)
      total=0
      

c     Compte le nombre de gal par bin de couleur
      do i=1,ngals_ada 
       if(dabs(ab_ada(i,1)).le.50.and.sab_ada(i,1).gt.0.)then
        if(dabs(ab_ada(i,4)).le.50.and.sab_ada(i,4).gt.0.)then
          x=ab_ada(i,1)-ab_ada(i,4)
          if(x.le.xmin)x=xmin
          if(x.ge.xmax)x=xmax
          xcmp(nint((x-xmin)/xstep))=xcmp(nint((x-xmin)/xstep))+1.d0
          total=total+1.d0
        endif
       endif
      enddo


c     fraction de gal par bin de couleur
      do i=0,nstep
        xcmp(i)=xcmp(i)/total
        write(112,*)xmin+dble(i)*xstep,xcmp(i)
      enddo


c     Attribut a chaque gal un weight
      do i=1,ngals_ada
       weight(i)=1.d0
       if(dabs(ab_ada(i,1)).le.50.and.sab_ada(i,1).gt.0.)then
        if(dabs(ab_ada(i,4)).le.50.and.sab_ada(i,4).gt.0.)then
          x=ab_ada(i,1)-ab_ada(i,4)
          if(x.le.xmin)x=xmin
          if(x.ge.xmax)x=xmax
          if(xcmp(nint((x-xmin)/xstep)).gt.0.d0)then
           weight(i)=1/xcmp(nint((x-xmin)/xstep))
           if(weight(i).gt.100)weight(i)=100.d0
           if(nint((x-xmin)/xstep).eq.0)write(*,*)"=0",i,weight(i)
          endif
          write(113,*)i,x,weight(i)
        endif
       endif        
      enddo




      end
