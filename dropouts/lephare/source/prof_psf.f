        program   prof
c
      integer*4 i,n 
      real*8    x,y,z,w,x2,dz,nt,beta,sigma
      real*8    io,g,rest,sum,sumo,intp,intpn,gtinf
      real*8    bessi0,gab,pi
      EXTERNAL  bessi0,gab
c
      pi=3.141592654
      open(1,file='prof_psf.dat')
      do nt = 1, 2
        if (nt.eq.1) then 
           n=1
           gtinf = 2.22618527   !G(infiny)=Gt(infini)= (2n)!/an**2n
        elseif (nt.eq.2) then
           n=4
           gtinf = 40320*pi/(7.67**(2*n)) !G(infiny)=Gt(infini)= (2n)!/an**2n 
        endif   
        dz= 0.1
        do w = 1, 301
          sigma=(w-1)*dz
          intp=0
          do y = 1, 301
            beta=(y-1)*dz
            sum=0
            if (sigma.gt.0) then
              do i = 1, 900
                z=(i-1)*dz
                sumo=sum
                g=gab(n,z)
                x=z*beta/sigma**2
                io=bessi0(x)
                x2=z*z/(2*sigma**2)
                rest=z*dexp(-1*x2)
                sum=sum+g*io*rest*dz
                if ( ((sum-sumo)/sum).le.1.e-7 ) goto 10
                if ( dlog10(io).gt.250) goto 20
              enddo
 10           sum=sum*1./sigma**2*dexp(-beta**2/2/sigma**2)
              goto 30
 20           sum=gab(n,beta)
 30           intp = intp + 2*pi*beta*sum*dz
            else
              sum=gab(n,beta)
              intp = intp + 2*pi*beta*sum*dz
            endif   
            sum=dlog10(sum)
            intpn=intp/gtinf
           write(1,'(i2,2(f6.2,1x),2(E12.6,1x))') n,beta,sigma,sum,intpn
           write(*,'(i2,2(f6.2,1x),2(E12.6,1x))') n,beta,sigma,sum,intpn
          enddo
        enddo
      enddo  
      close(1)
      stop
      end
c
ccccccccccccccccccccccccccccccccccc
        REAL*8  FUNCTION gab(n,b) 
        implicit none 
        integer*4 n
        real*8    a,b 
c
        if (n.eq.1) a=1.68
        if (n.eq.4) a=7.67
        gab=dexp(-1.*a*(b**(1./n)))
c
        return
        end       
c
ccccccccccccccccccccccccccccccccccccccccccccc
      REAL*8 FUNCTION bessi0(x)
      REAL*8 x
      REAL*8 ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     *1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     *0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,
     *-0.1647633d-1,0.392377d-2/
      if (abs(x).lt.3.75) then
        y=(x/3.75)**2
        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax=dabs(x)
        y=3.75/ax
        bessi0=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *(q7+y*(q8+y*q9))))))))
      endif
      return
      END
