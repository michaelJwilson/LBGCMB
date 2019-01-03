c        15/04/03
c        Get two tables in the command line
c        For the same index, we sort tab1(i)<tab2(i)
c        The two table are double
cccccccccccccccccccccccccccccccccc
         subroutine getOpt2Tab(nbBin,tab1,tab2,lenght,
     .                         test,paravr,param)


c        nbBin* : number of bin deduced from command line (normal case test/2)
c        tab1*  : first output table (double)
c        tab2*  : second output table (double)
c        lenght : declaration lenght of both tables
c        test   : number of integer in command line
c        paravi : tab with integer read in command line
c        param  : name of option line command 
 
         implicit none

         INCLUDE 'out_unit.decl' 

         integer*4  i,lnblnk,nbBin,lenght,test
         real*8 tab1(lenght),tab2(lenght),paravr(500)
         character*4096 param
         
         
c        If you want to give a name to the output screen file
         if(UO.eq.30)
     .      open(UO,file='/tmp/screenGet2Tab.dat',status='unknown') 


         if(test.le.1)then
           write(UO,*)param(1:lnblnk(param))," not defined" 
           do i=1,lenght
            tab1(i)=0
            tab2(i)=100
           enddo
         elseif(test.eq.2)then
            nbBin=1
            if(paravr(1).lt.paravr(2))then
              do i=1,lenght
               tab1(i)=paravr(1)
               tab2(i)=paravr(2)
              enddo
             else
              do i=1,lenght
                tab1(i)=paravr(2)
                tab2(i)=paravr(1)
              enddo
             endif
         else
            nbBin=nint(dble(test)/2.d0)
            do i=1,test-1,2
              if(paravr(i).lt.paravr(i+1))then
               tab1(nint(dble(i)/2.))=paravr(i)
               tab2(nint(dble(i)/2.))=paravr(i+1)
              else
               tab1(nint(dble(i)/2.))=paravr(i+1)
               tab2(nint(dble(i)/2.))=paravr(i)
              endif
            enddo
         endif

         return


         end 
