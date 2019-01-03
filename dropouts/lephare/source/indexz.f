c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      INTEGER*4 FUNCTION  indexz(z,zstep,dz)
      implicit none
      real*8   z,zstep,dz
c
      if (z.le.6) then
         indexz = idnint(z/zstep) + 1
      else
         indexz = idnint(6./zstep) + 1 + idnint((z-6.)/dz)
      endif   
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

