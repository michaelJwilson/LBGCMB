c
cccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine kInterp(zs,index,imagm,zlib,modlib,inputLib,outTab)
c
c     Author: Olivier Ilbert (20/03/03)
c     Goal:   Linear interpolation of k correction 
c
      implicit none
      integer*4 nbf,inlib
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_lib.decl'
      integer*4 index,sens,imagm,k,modlib(inlib)
      real*8 zlib(inlib),inputLib(nbf,inlib),outTab(nbf),zs,a
c
c     sens permet de se positionner dans le bin encadrant zs     
      sens=idnint(dsign(1.d0,(zs-zlib(index))))     
c
c     linear factor for interpolation
      a=dabs((zlib(index)-zs)/(zlib(index)-zlib(index+sens)))
c
c     Interpolation for each filter the kcorr or app. mag.       
      do k=1,imagm
         if(modlib(index).ne.modlib(index+sens))      a=0.d0
         if( inputLib(k,index+sens).gt.90)            a=0.d0
         outTab(k) = a*inputLib(k,index+sens)+(1-a)*inputLib(k,index)
      enddo
c
      return
      end
