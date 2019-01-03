ccccccccccccccccccccccccccccccccccccc
        real*8 function dmet(z,h0,om0,l0)
ccccccccccccccccccccccccccccccccccccc
c      Compute the metric distance dmet in Mpc
c      dlum = dmet*(1+z) 
c      dang = dmet/(1+z) = dlum/(1+z)^2
c
        implicit none
        integer*4 i 
        real*8    om0,h0,z,c,l0,omt,ao,sum,dz,zi,Ez
c       
        dmet=0.
        c = 300000.
        omt = om0+l0
        if (l0.eq.0) then
c          ao = c/(h0*sqrt(ABS(1-omt)))
c  in fact we use x = ao * x(z) with x(z) from eq 8 of 
c    Moscardini et al.  So we don't need to compute ao   
          if (om0.gt.0) then
            ao =1.
            dmet = om0*z-(om0-2)*(1-dsqrt(1+om0*z))
            dmet= 2*c/(ao*h0*om0*om0*(1+z))*dmet
          else
            ao =1.
            dmet= c*z*(1+z/2)/(h0*(1+z))
          endif 
        elseif (om0.le.1.and.l0.ne.0) then
          ao = 1.
          sum = 0. 
          dz = z/50.
          do i = 1,51
            zi = (i-0.5)*dz
            if (zi.gt.z) goto 1
            Ez= dsqrt(om0*(1+zi)**3+(1-om0-l0)*(1+zi)**2+l0)
            sum = sum + dz/Ez
          enddo        
 1        dmet =  c/(h0*ao) * sum
        endif         
        return
        end
ccccccccccccccccccccccccccccccccccccc



