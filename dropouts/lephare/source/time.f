c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Compute cosmological time from z=infinty  to z
c     as a function of cosmology.  Age given in year !!
c     Note : use lambda0 non zero if Omega_o+lambda_o=1
	real*8 function timy(z,h0,om0,l0)
	implicit none
	REAL*8  hy,h0,l0,z,om0,val
c
c   Gives time in yr
	timy=0.
        hy=h0*1.0224e-12
        if (dabs(om0-1).le.1.e-6 .and. l0.eq.0) then
           timy=2*(1+z)**(-1.5)/(3*hy)
        elseif(om0.eq.0.and.l0.eq.0) then
           timy=1./(hy*(1+z))
        elseif (om0.lt.1.and.om0.gt.0 .and. l0.eq.0) then
	   val=(om0*z-om0+2)/(om0*(1+z))
           timy=2.*dsqrt((1-om0)*(om0*z+1))/(om0*(1+z))
           timy=timy - dlog(val+dsqrt(val*val-1))
           timy=timy*om0/(2.*hy*(1-om0)**(1.5))
        elseif (om0.gt.1 .and. l0.eq.0) then
           timy=dacos((om0*z-om0+2)/(om0*(1+z)))
           timy=timy-2*dsqrt((om0-1)*(om0*z+1))/(om0*(1+z))
           timy=timy*om0/(2*hy*(om0-1)**1.5)       
        elseif (om0.lt.1 .and. DABS(om0+l0-1).lt.1.e-5 ) then
	   val =  dsqrt(1-om0)/(dsqrt(om0)*(1+z)**1.5)
	   timy= dlog(val+dsqrt(val*val+1))
           timy= timy*2./(3.*hy*dsqrt(1-om0))
        endif     
	return
	END
c
