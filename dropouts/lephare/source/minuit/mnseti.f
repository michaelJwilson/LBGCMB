*
* $Id: mnseti.F,v 1.1.1.1 1996/03/07 14:31:31 mclareni Exp $
*
* $Log: mnseti.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
c#include "minuit/pilot.h"
      SUBROUTINE MNSETI(TIT)
      include 'd506dp.inc'
CC       Called by user to set or change title of current task.
CC
      include 'd506cm.inc'
      CHARACTER*(*) TIT
      CTITL = TIT
      RETURN
      END
