c     last modif 21/11/02
c     description : Return character value for floating conversion 
ccccccccccccccccccccccccccccccc
      subroutine check_float(string,str_ch)
cccccccccccccccccccccccccccccc
      implicit none
      integer*4     i,lnblnk,imin,idot
      character*4096 string,var,str_ch
c
      str_ch= string(1:lnblnk(string))
c  check if a "." is detected to avoid problem in float conversion
      imin=0
      idot=0
      do i = 1, lnblnk(string)+1
               var=string(i:i)
               if (var(1:1).eq.'.') then
                  idot=1
               endif   
               if (var(1:1).eq.'e' .or. var(1:1).eq.'E'
     >        .or. var(1:1).eq.'d' .or. var(1:1).eq.'D') imin=i 
            enddo
            if (idot.eq.0 .and.imin.eq.0) then
                str_ch= string(1:lnblnk(string)) // '.'
            elseif (idot.eq.0 .and.imin.ne.0) then 
               str_ch =  string(1:imin-1) 
     >             // '.' // string(imin:lnblnk(string))
            endif   
c 
      RETURN
      END
