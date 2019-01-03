C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        real*8 function zeta(t,h0,om0,l0)
        implicit none
c
        real*8 timy,zmin,zmax,tmax,tmin,t0
        real*8 zmed,tmed,dt,t,h0,l0,om0,z0
c 
        external  timy
c
	zeta = 0.
        z0=0.
        if (t.le.1.e-7) return 
        zmin = 0.
        zmax = 1000.
        t0 = timy(z0,h0,om0,l0)
        tmin = t0-timy(zmin,h0,om0,l0)
        tmax = t0-timy(zmax,h0,om0,l0)
        dt = tmax - tmin
        do while (dt .ge. 0.00001*1.E9)
          zmed = (zmin + zmax) / 2.
          tmed = t0-timy(zmed,h0,om0,l0)
          if (t.le.tmax .and. t.gt.tmed) then
             zmax = zmax
             zmin = zmed
          else if (t.le.tmed .and. t.ge.tmin) then          
             zmax = zmed
             zmin = zmin
          endif
          tmax = t0-timy(zmax,h0,om0,l0)
          tmin = t0-timy(zmin,h0,om0,l0)
          dt = tmax - tmin
        end do
        zmed = (zmin + zmax) / 2.
        zeta = zmed       
	return
	END
