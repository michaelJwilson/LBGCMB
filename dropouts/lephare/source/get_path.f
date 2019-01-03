c      last modif 30/05/2003
c      Author : S. ARNOUTS 
c      description : Extract environment pathname from character string
c      if not found returns test=0
ccccccccccccccccccccccccccccccccccccccc 
      subroutine get_path(string)
      implicit none 
c     
      INCLUDE 'out_unit.decl' 
      integer*4     i,imin,imax,ilen,lnblnk
      character*4096  string,path
      character*4096  var
c
c      If you want to give a name to the output screen file
      if(UO.eq.30)
     .      open(UO,file='/tmp/screenGatPath.dat',status='unknown') 
c
c     search the first non blank character and the last character
c      write(6,*) 'in :',string(1:lnblnk(string)) 
      if ( string(1:1) .eq. '$' ) then 
         imin=1
         imax=0
         do i = 2, 4096
           var = string(i:i) 
           if ( var(1:1).eq.'/' .or. var(1:1).eq.'#'
     >         .or. var(1:1).eq.char(9) .or. var(1:1).eq.' ') then
              imax=i-1
              goto 1
           endif
         enddo
 1       if (imax.le.imin) then
            write(UO,*) ' ERROR : Failed to recognize the pattern: ',
     >                  string(1:lnblnk(string)),' --> STOP'
            stop
         else
            var=string(imin+1:imax)
c            write(UO,*) var(1:lnblnk(var)),string(imin+1:imax)
            call getenv(var(1:lnblnk(var)),path)
            
            if (lnblnk(path) .eq. 0) then
      write(UO,*) ' ERROR: $',path(1:lnblnk(path)),' not found --> STOP'
              stop
            else
               ilen = lnblnk(string) - imax
               string = path(1:lnblnk(path)) // 
     >                  string(imax+1:lnblnk(string))
c               write(UO,*) 'out :',string(1:lnblnk(string)) 
            endif
            
         endif   
      endif   
c      write(UO,*) ' -> return: ',string(1:lnblnk(string)) 
c
      return
      end
