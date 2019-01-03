c
cccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine z_vmax(index,dmag,zbd,modlib,extlib,ebvlib,
     > agelib,klib,zlib,dmlib,lmax,zmaxlib)
c
c     Author: S. Arnouts 
c     Goal:   extract Zmax for a given SED knowing the mag_lim 
c      Diff=Mag_lim-Mag_obs=(DM+kcor)(zmax)-(DM+kcor)(zobs)
c      -> Input      Dmag= Diff+(DM+kcor)(zobs) 
c      -> search for Dmag=(DM+kcor)@ z=zmax        
      implicit none
      integer*4 nbf,inlib
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_lib.decl'
      integer*4 index,lmax,zbd,i
      integer*4 modlib(inlib),extlib(inlib),modref,extref
      real*8    ebvlib(inlib),agelib(inlib),dmlib(inlib),zlib(inlib)
      real*8    klib(nbf,inlib),ebvref,ageref
      real*8    dmag,diff,zmaxlib,diffp,zp
c
c   assumes the model,Extinction, age at index 
      modref=modlib(index)
      extref=extlib(index)
      ebvref=ebvlib(index)
      ageref=agelib(index)
c
c      write(*,*) modref,zbd
c
      zmaxlib=99.99
      diffp=0
      zp=zlib(index)
      do i = index,lmax
         diff = dmag - (dmlib(i)+klib(zbd,i))
         
c         write(*,'(3(I8,2x),4(E12.6,2x),2(f6.3,2x))')
c     >       i,modlib(i),extlib(i),ebvlib(i),agelib(i),
c     >       dmlib(i),klib(zbd,i),zlib(i),diff

         if (i.eq.index .and. diff.lt.0) then 
             zmaxlib=zlib(i)
             goto 10 
         endif
 
         if (diff.gt.0) then 
            diffp=diff
            zp = zlib(i)
         endif
         if (i.gt.index .and. diff.lt.0) then 
             zmaxlib=zlib(i)
             goto 10 
         endif
      
         if ( modlib(i) .ne. modref .or. 
     >        extlib(i) .ne. extref .or. 
     >        ebvlib(i) .ne. ebvref .or. 
     >        agelib(i) .ne. ageref      ) then 
             zmaxlib=zlib(i-1)
             goto 10 
         endif 
      enddo 
 10   continue
c      write(*,*) 'Zmax =',zbd,diffp,zp,diff,zmaxlib
      if (diff.lt.0 .and. diffp.gt.0) then 
          zmaxlib= (diff*zp-diffp*zmaxlib)/(diff-diffp)         
      endif
c      write(*,*) 'Zmax =',zbd,zmaxlib  
c
      return
      end
