c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Return zmin,zmax for chi2+Dchi2 (68%,90%,99%)
c
      SUBROUTINE DCHI2(chi2,bmax,chib,dchi,zinf,zsup)
c
      implicit none 
c      
      INTEGER*4 bmax,i,k
      REAL*8    chi2(2,bmax),chib,dchi,chi,slope
      REAL*8    zi(1000),zs(1000),zinf,zsup
c
c   Initialisation
      chi=chib+dchi
      k = 0 
      if (chi2(2,1).le.chi) then
         k    = k + 1
         zi(k) = chi2(1,1)
      ENDIF   
      do i = 1,bmax-1
         if (chi2(2,i).gt.chi .and. chi2(2,i+1).le.chi) then
            k     = k + 1
c            zi(k) = chi2(1,i)
            slope= (chi2(2,i+1)-chi2(2,i)) / (chi2(1,i+1)-chi2(1,i))
            zi(k) = chi2(1,i)+ (chi-chi2(2,i))/slope

         endif
         if (chi2(2,i).lt.chi .and. chi2(2,i+1).ge.chi) then
c            zs(k) = chi2(1,i+1)
            slope= (chi2(2,i+1)-chi2(2,i)) / (chi2(1,i+1)-chi2(1,i))
            zs(k) = chi2(1,i)+ (chi-chi2(2,i))/slope           
         endif
      enddo   
      if (chi2(2,bmax).le.chi) then
         zs(k) = chi2(1,bmax)
      ENDIF        
c
      if (k.eq.1) then
         zinf=zi(1)
         zsup=zs(1)
      elseif (k.gt.1) then 
         zinf=zi(1)
         zsup=zs(k)
      elseif (k.eq.0) then
         zinf=0
         zsup=0
      endif   
c   
      return
      END
c
