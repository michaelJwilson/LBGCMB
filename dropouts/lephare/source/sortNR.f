cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       sort using a numerical recipes subroutine

        subroutine sortNR(tab,nb1,nb2,index)

c       tab*    : input table to sort
c       nb1     : lenght of tab to sort
c       nb2     : lenght of the input array
c       index   : index of the values in tab sorted by increasing values

        implicit none

        integer*4 nb1,nb2
        real*8 rra,tab(nb2),ra2(nb2)
        integer*4 i,j,l,ir,iid,index(nb2)



c       Copy the input array, I don't want to overwrite it        
        do i=1,nb2
          ra2(i)=tab(i)
          index(i)=i
        enddo

        l=nb1/2+1
        ir=nb1
   10   continue
        if(l .gt. 1)then
         l=l-1
         rra=ra2(l)
         iid=index(l)
        else
         rra=ra2(ir)
         iid=index(ir)
         ra2(ir)=ra2(1)
         index(ir)=index(1)
         ir=ir-1
         if(ir.eq.1)then
          ra2(1)=rra
          index(1)=iid
          return
         endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
         if(j .lt. ir)then
          if(ra2(j) .lt. ra2(j+1))  j=j+1
         endif
         if(rra .lt. ra2(j))then
          ra2(i)=ra2(j)
          index(i)=index(j)
          i=j
          j=j+j
         else
          j=ir+1
         endif
         goto 20
        endif
        ra2(i)=rra
        index(i)=iid
        goto 10
       end

        
        

