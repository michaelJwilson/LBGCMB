      subroutine zmax(zs,ab,imagm,zlib,modlib,klib,vmax)
c
c     Goal:   Maximum Z knowing the maglim in a given band 
c     mlim = mag +(kcor(zmax)-kcor(z0)) + (DM(zmax)-DM(z0)) 
      implicit none
      integer*4 nbf,inlib
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_lib.decl'
      integer*4 index,sens,imagm,k,modlib(inlib)
      real*8 zlib(inlib),kib(nbf,inlib),outTab(nbf),zs,a
      real*8 ab(nbf)
c




c
      return
c
