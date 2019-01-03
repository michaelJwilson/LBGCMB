cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE sort_col(flux,nf,max,ncol,color,c_index)
      implicit none 
      integer*4  max,nf,inlib,nbf
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_lib.decl'
      integer*4  indx(inlib),c_index(nbf,inlib)
      integer*4  i,k,ncol
      real*8     flux(nbf,inlib)
      real*8     color(nbf,inlib),col_k(inlib)
c
c      write(6,*) 'sorting the colors ... '
      do k = 1,ncol
         do i = 1,max
            if (flux(k+1,i).gt.0) then
               col_k(i)= flux(k,i)/flux(k+1,i)
            else
               col_k(i) = 0.
            endif   
         enddo
         call  indexx(max,col_k,indx)
         do i = 1,max
           color(k,i)= col_k(indx(i))
           c_index(k,i)=indx(i)
         enddo  
      enddo
c      do i = 1, max
c           write(*,*) i,(color(k,i),c_index(k,i),k=1,ncol)
c      enddo     
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
