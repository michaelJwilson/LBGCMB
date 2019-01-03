c        Get one tables from command line
c        The table contain interger value
c        15/04/03
         subroutine getOpt1Tab(nbBin,tab,lenght,test,paravi,param)

c        nbBin* : number of bin deduced from command line
c        tab*   : table with integer of command line
c        lenght : declaration lenght of tab
c        test   : number of integer in command line
c        paravi : tab with integer read in command line
c        param  : name of option line command 
 
         implicit none

         INCLUDE 'out_unit.decl' 

         integer*4  i,lnblnk,nbBin,lenght,test
         integer*4 tab(lenght),paravi(500)
         character*4096 param
         
c        If you want to give a name to the output screen file
         if(UO.eq.30)
     .      open(UO,file='/tmp/screenGet1Tab.dat',status='unknown') 


         if(test.eq.0)then
c          If no parameter read, put 1 in tab
           write(UO,*)param(1:lnblnk(param))," not defined"
            nbBin=1
           do i=1,lenght
             tab(i)=1
           enddo
         elseif(test.eq.1)then
c           If only one parameter read, put this same value in all bin
            nbBin=1
            do i=1,lenght
             tab(i)=paravi(1)
            enddo
         else
            nbBin=test
            do i=1,test
               tab(i)=paravi(i)
            enddo
         endif

         return


         end 
