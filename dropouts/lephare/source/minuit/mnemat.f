*
* $Id: mnemat.F,v 1.1.1.1 1996/03/07 14:31:29 mclareni Exp $
*
* $Log: mnemat.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
c#include "minuit/pilot.h"
      SUBROUTINE MNEMAT(EMAT,NDIM)
      include 'd506dp.inc'
      DIMENSION EMAT(NDIM,NDIM)
CC        Calculates the external error matrix from the internal
CC        to be called by user, who must dimension EMAT at (NDIM,NDIM)
      include 'd506cm.inc'


cccccccccccccccccccc
c       J'ajoute
      common/ErrorMatrixSTY/EMATCOMSTY(2,2)
      common/ErrorMatrixFIT/EMATCOMFIT(3,3)
cccccccccccccccccccc

      IF (ISW(2) .LT. 1)  RETURN
      IF (ISW(5) .GE. 2)  WRITE (ISYSWR,'(/A,I4,A,I3,A,G10.3)')
     +    ' EXTERNAL ERROR MATRIX.    NDIM=',NDIM,'    NPAR=',NPAR,
     +    '    ERR DEF=',UP
C                    size of matrix to be printed
      NPARD = NPAR
      IF (NDIM .LT. NPAR)  THEN
        NPARD = NDIM
        IF (ISW(5) .GE. 0) WRITE (ISYSWR,'(A,A)') ' USER-DIMENSIONED ',
     +      ' ARRAY EMAT NOT BIG ENOUGH. REDUCED MATRIX CALCULATED.'
      ENDIF
C                 NPERLN is the number of elements that fit on one line
      NPERLN = (NPAGWD-5)/10
      NPERLN = MIN(NPERLN,13)
      IF (ISW(5).GE. 1 .AND. NPARD.GT.NPERLN)  WRITE (ISYSWR,'(A)')
     +     ' ELEMENTS ABOVE DIAGONAL ARE NOT PRINTED.'
C                 I counts the rows of the matrix
      DO 110 I= 1, NPARD
         CALL MNDXDI(X(I),I,DXDI)
         KGA = I*(I-1)/2
         DO 100 J= 1, I
            CALL MNDXDI(X(J),J,DXDJ)
            KGB = KGA + J
            EMAT(I,J) = DXDI * VHMAT(KGB) * DXDJ * UP
            EMAT(J,I) = EMAT(I,J)
  100    CONTINUE
  110 CONTINUE
C                    IZ is number of columns to be printed in row I
      IF (ISW(5) .GE. 2)  THEN
      DO 160 I= 1, NPARD
         IZ = NPARD
         IF (NPARD .GE. NPERLN)  IZ = I
         DO 150 K= 1, IZ, NPERLN
           K2 = K + NPERLN - 1
           IF (K2 .GT. IZ)  K2=IZ
           WRITE (ISYSWR,'(1X,13E10.3)')  (EMAT(I,KK),KK=K,K2)
  150    CONTINUE
  160 CONTINUE
      ENDIF
cccccccccccccccccccc
c     AJOUTE     
cccccccccccccccccccc
c     AJOUTE     
      if(NPARD.eq.3)then
       do i=1,3
        do j=1,3
          EMATCOMFIT(i,j)=EMAT(i,j)
        enddo
       enddo
      endif
      if(NPARD.eq.2)then
       do i=1,3
        do j=1,3
          EMATCOMSTY(i,j)=EMAT(i,j)
        enddo
       enddo
      endif
cccccccccccccccccccc
      RETURN
      END
