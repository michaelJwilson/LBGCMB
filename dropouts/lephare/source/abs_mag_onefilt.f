cccccccccccccccccccccccccccccccccccccccc
c  Compute absolute magnitude
c
      subroutine absMag1F(
     >    method,Id,fabs,bapp,zstep,dz,moddist,
     >    filter_obs,filter_obs4,bused,z,dm,ab,magm,sab,kcor,maglib0,
     >    min_kcolor,goodfilter,minkcolor,pdz_mabsz)
c
      implicit none
c
      integer*4 k,l,nbf,chisize
c
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_zchi.decl'
c
      integer*4  method(nbf),methodold,fabs,bapp(nbf),Id
      real*8     errLim,z,moddist,dm,zstep,dz
      integer*4  filter_obs(chisize,nbf,nbf)
      integer*4  filter_obs4(chisize,nbf,nbf)
      integer*4  goodfilter(nbf),Fobs,bused(nbf,nbf),findfilter
      real*8     ab(nbf),sab(nbf),magm(nbf),kcor(nbf)
c      real*8     sabm(nbf),zp(nbf)
c      real*8     cab(nbf)
c      integer*4  inlib
      real*8     min_kcolor(chisize,nbf,nbf),minkcolor(nbf),maglib0(nbf)
      real*8     color,pdz_mabsz,kcolor
c      real*8     funz,h0,om0,l0
c
      integer*4    indexz
      real*8       d_mod
      external      funz,d_mod,indexz
c
c      moddist=funz(z,h0,om0,l0) !Compute DM at accurate z
c
      do k=fabs,fabs !Compute M in all filter

        methodold=-1
        errLim=0.6
 6      goodfilter(k)=1
        minkcolor(k)=9999
        findfilter=0
c
c   method = 0 : used observed mag and convert

        if (method(k) .eq. 0) then 
           goodfilter(k)=k
        endif   
c
c   method 1 : select filter closed to the rest-frame reference filter
        if(method(k).eq.1 .or. method(k).eq.4)then
          do l=1,nbf
           if (method(k) .eq. 1) then 
              Fobs=filter_obs(indexz(z,zstep,dz),k,l) ! Select best filter
           else
              Fobs=filter_obs4(indexz(z,zstep,dz),k,l) ! Select best filter
           endif
           if (bused(k,Fobs).eq.1)then
            if(sab(Fobs).ge.0.and.sab(Fobs).le.errLim)then
             if(ab(Fobs).lt.90)then
              if(kcor(Fobs).lt.90)then
                goodfilter(k)=Fobs
                findfilter=findfilter+1
                minkcolor(k)=min_kcolor(indexz(z,zstep,dz),k,l)
                goto 7
              endif
             endif
            endif
           endif
          enddo
c
c          If unable with constrain in error, try without errlim
          if(findfilter.eq.0.and.errLim.lt.999.9)then  
c           write(6,*)"Enable to find adapted filter for ",Id
c           write(6,*)"try without err lim              "
           errLim=99999999999.99
           goto 6
          endif

c          write(6,*)"Can t use method 1 for : ",Id," Use method 3   "
          methodold=method(k)
          method(k)=-1         
c
        endif

c
c   Method 2 : Fix bapp as observed-frame
        if(method(k).eq.2)then
         if(bused(k,bapp(k)).eq.1.and.ab(bapp(k)).lt.90
     .   .and.sab(bapp(k)).ge.0.and.kcor(bapp(k)).lt.90)then
           goodfilter(k)=bapp(k)
           goto 7
         else
c          write(6,*)"Can t use method 2 for : ",Id," Use method 3   "
           methodold=method(k)
           method(k)=-1
         endif
        endif


 7      continue

c       Compute absolute magnitude
        if( method(k).eq.0 .or. method(k).eq.1
     >     .or. method(k).eq.2  .or. method(k).eq.4.)then    

c         Method 1 and 2 : Use observed apparent magnitude
c         For method 1 : best filter
c         For method 2 : Imposed filter
          color=maglib0(k)-maglib0(goodfilter(k)) !Rest-frame color
          kcolor=color-kcor(goodfilter(k))        !rf color-kcor
          pdz_mabsz=ab(goodfilter(k))+kcolor-moddist !AbsM from observed mag
          
        elseif (method(k).eq.3 .or. method(k).eq.-1) then

c         Absolute magnitude with method 3.
c         Use integration of the flux from the SED * scaling for observed flux
c          abs_mag(k)=magm(k)-moddist-kcor(k)
c          abs_mag(k)=abs_mag(k)-2.5*dlog10(dm)
           pdz_mabsz = maglib0(k) -2.5*dlog10(dm) 
          if (method(k) .eq.  3)  goodfilter(k) = 0
          if (method(k) .eq. -1)  then 
             method(k)=methodold
             goodfilter(k) = -1
          endif   
        endif
c
       enddo
c
       return
       end
