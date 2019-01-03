c   Distance modulus given in Mpc 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 FUNCTION  FUNZ(Z,H0,om0,l0)
      IMPLICIT NONE 
      REAL*8  Z,H0,om0,l0,dmet
c
      external  dmet
c
c  Dlum = Dmet*(1+z)
      IF ( Z.le.1.D-5 ) then
          funz = 0 
      ELSE
         funz=5.*dlog10((1+Z)*dmet(Z,H0,om0,l0)) + 25
      ENDIF
c
      RETURN
      END
