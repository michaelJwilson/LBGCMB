c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine lbref(lref,iw,lini,fini,nini,fout)
      implicit none
c
      integer*4  i,j,nini,iw,wmax
c      parameter (wmax=8000)
      INCLUDE 'dim_wave.decl'      
      real*8     lref(wmax),lini(wmax),fini(wmax),fout(wmax)
      real*8     lam,dlamb
c
c  convertion in the same lbda as wave from model
c       write(*,*) iw
       do i = 1,iw
          fout(i) = 0.
          lam = lref(i)
          if (lam.lt.lini(1).or.lam.gt.lini(nini)) then
             fout(i) = 0.
          else
             do j = 1,nini-1
                dlamb = lam - lini(j)
                if (dabs(dlamb).le.0.005) then
                   fout(i) = fini(j)
                   goto 10
                elseif (lam.gt.lini(j).and.lam.lt.lini(j+1)) then 
                  fout(i) = (fini(j) + fini(j+1)) / 2.
                   goto 10
                endif
             enddo
             dlamb = lam - lini(nini)
             if (dabs(dlamb).le.0.005) then
                fout(i) = fini(nini)
             endif
             goto 10                
 10           continue
          endif
        enddo                    
c
      return
      end
c
 
