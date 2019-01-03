c
      subroutine err_option(var,num)
      implicit none 
      integer*4     num,lnblnk
      character*4096 var
c
      if (num.eq.1) then
        write(*,*) 'Main Parameter "',var(1:lnblnk(var)),
     $            '"  not found -> STOP'
        STOP
      elseif (num.eq.2) then
        write(*,*) 'Option "',var(1:lnblnk(var)),
     $            '" not defined'
      elseif (num.eq.3) then 
        write(*,*) 'Option "',var(1:lnblnk(var)),
     $            '" not the right number of parameters  -> STOP' 
        STOP
      elseif (num.eq.4) then 
        write(*,*) 'Option "',var(1:lnblnk(var)),
     $            '" not the right number of parameters' 
      endif
c
      RETURN
      END
