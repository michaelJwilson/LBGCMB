*
* $Id: d506dp.inc,v 1.1.1.1 1996/03/07 14:31:32 mclareni Exp $
*
* $Log: d506dp.inc,v $
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
c#ifndef CERNLIB_MINUIT_D506DP_INC
c#define CERNLIB_MINUIT_D506DP_INC
*
*
* d506dp.inc
*
c#if defined(CERNLIB_DOUBLE)
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c#endif
c#if !defined(CERNLIB_DOUBLE)
C ************ SINGLE PRECISION VERSION *************
c#endif

c#endif
