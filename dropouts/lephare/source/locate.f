cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE locate(xx,n,x,j)
      INTEGER*4 j,n
      REAL*8    x,xx(n)
      INTEGER*4 jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END
c
