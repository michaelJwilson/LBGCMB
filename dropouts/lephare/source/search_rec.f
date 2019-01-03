c   
      SUBROUTINE sear_rec(fl,efl,nf,band,sig,ncol,selec,
     >                colt,indt,tmax,rec_read,nrec)
c
      implicit none 
c
      integer*4  nf,inlib,nbf,ncol,nrec,tmax,nused
      INCLUDE 'dim_filt.decl'
      INCLUDE 'dim_lib.decl'
c
      integer*4 band(nbf),indt(nbf,inlib)
      integer*4 rec_read(inlib),irec(nbf),rec_col(nbf,inlib)
      real*8 fl(nbf),efl(nbf)
      real*8 sig,colt(nbf,inlib),colort(inlib)
      real*8 colmin,colmax
      real*8 color,errcol
      character*4096 selec
c
      integer*4  i,k,imin,imax
c 
c
      color=0      
      imin=0 
      imax=0
c Check if colors used 
      nused=0
      do k = 1, ncol
         if((band(k)*band(k+1)).eq.1) nused=nused+1 
      enddo
      if (nused.eq.0) then
         do i = 1 , tmax
          rec_read(i)=indt(1,i)        
         enddo
         nrec=tmax
         return
      endif
c  Otherwise search color interval   
      do k = 1, ncol 
c  Check if upper-limits (err_fl=-1) -> err_fl=fl_lim 
         if (efl(k).eq.-1)   efl(k)  = DABS(fl(k))   ! In any case, must be 
         if (efl(k+1).eq.-1) efl(k+1)= DABS(fl(k+1)) ! already positive !!
         do i = 1,tmax
            rec_col(k,i) = 0
            colort(i)    = colt(k,i)
         enddo
         if((band(k)*band(k+1)).eq.0) then      ! At least 1 band not used 
            goto 1
         elseif ( (band(k)*band(k+1)).eq.1 .AND.
     >             fl(k+1).gt.0 .AND. fl(k).gt.0 ) then
             color = fl(k)/fl(k+1)
         elseif ( (band(k)*band(k+1)).eq.1 .AND.
     >       ( fl(k+1).le.0 .OR. fl(k).le.0 ) )  then
             color = 0.
         endif   
c      write(*,*) 'search record 1 ',k
         if (fl(k+1).gt.0) then
           errcol = (efl(k)*fl(k+1))**2
           errcol = errcol + (efl(k+1)*fl(k))**2
           errcol = dsqrt(errcol)/(fl(k+1)**2)
         else
           errcol = (efl(k)*fl(k+1))**2
           errcol = errcol +(efl(k+1)*fl(k))**2
           errcol = dsqrt(errcol)/(efl(k+1)**2) 
         endif  
         colmin = color-sig*errcol
         colmax = color+sig*errcol
         irec(k)= 0
        call locate(colort,tmax,colmin,imin)
c        call locate(colort,inlib,colmin,imin)
        imin=imin+1
        call locate(colort,tmax,colmax,imax)
c        call locate(colort,inlib,colmax,imax)
        do i = imin,imax
           irec(k) = irec(k)+1
           rec_col(k,indt(k,i)) = indt(k,i)
        enddo           
 1      continue
c         write(*,*) 'search record  ',k,colmin,colmax,imin,imax
      enddo    
c
c   SEARCH FOR COMMON RECORDS
      nrec=0
      if (selec(1:3) .eq. 'AND') then   
        do i = 1, tmax
          do k = 1,ncol
            if((band(k)*band(k+1)).eq.0) goto 20
            if (rec_col(k,i).eq.0) goto 21
 20         continue
          enddo
          nrec = nrec+1
          rec_read(nrec) = rec_col(1,i)        
 21       continue
        enddo  
      elseif (selec(1:2) .eq. 'OR') then 
        do i = 1, tmax
          do k = 1,ncol
             if (rec_col(k,i).ne.0) then
                nrec = nrec+1
                rec_read(nrec) = rec_col(k,i)     
                goto 22
             endif   
          enddo
 22       continue
        enddo           
      endif
      if (nrec.eq.0) then
         do i = 1 , tmax
          rec_read(i) = indt(1,i)        
         enddo
         nrec = tmax
      endif   
c
      return
      end
c



