      subroutine prep_out(cattyp,cat,config,outpara,imagm,
     >  wpara,iwpara,npara,nppara)
c
      implicit none  
      integer*4       test,conti,i,j,k,imagm,npdz
      integer*4       iwpara,imin,imax
      integer*4       npara(500),nppara(50)
      integer*4       paravi(500)
      real*8          val1,paravr(500)
      character*4096  cat,file,config,str,cattyp,param
      character*4096  paravc(500),outpara,wpara(500)
      character*4096  str1,str2,str3,ich
      INCLUDE 'out_unit.decl'
c
c  read cat if type == LONG for length of string_input
      conti=0
      if (cattyp(1:4).eq.'LONG') then
        file = cat 
        open(1,file=file,status='unknown',err=56)      
        do while (.true.)              
           read(1,'(A)',end=18) str 
           call val_string(str,paravc,test)
           if (paravc(1)(1:1) .ne. '#') then  
              j = 2*imagm+3
              if (test.gt.j) then
                conti= test-j
              else
                conti=0
              endif
              goto 18 
          endif
        enddo   
 18     close(1)
      endif
      npdz=0
      do i = 1,50
         nppara(i)=0
      enddo
ccccccccccccccccc
      file = outpara(1:lnblnk(outpara)) 
      open(1,file=file,status='unknown',err=56)
      iwpara = 0 
      do while (.true.)
         read(1,'(A)',end=13) str 
         if (str(1:1) .ne. '#' .AND. str(1:1) .ne. ' ' 
     >  .AND. str(1:1) .ne. char(9) .AND. str(1:1) .ne. char(13)) then
           call val_string(str,paravc,test)
            iwpara = iwpara + 1
            wpara(iwpara) = paravc(1)
c
          if (paravc(1)(lnblnk(paravc(1))-1:lnblnk(paravc(1))).eq.'()'
     >        .and. paravc(1)(1:lnblnk(paravc(1))).ne.'PDZ()'
     >        .and. paravc(1)(1:lnblnk(paravc(1))).ne.'Z_MAX()'  )then     
 
              if (iwpara.gt.1) then 
                 npara(iwpara)=imagm+npara(iwpara-1)
              else
                 npara(iwpara)=imagm
              endif
c
          elseif (paravc(1)(1:lnblnk(paravc(1))).eq.'PDZ()')then     

              param='-PROB_INTZ'
              call getf_option(param,config,500,paravr,test)
              if (test.ge.1 .and. MOD(test,2).eq.0) then
                 if (iwpara.gt.1) then 
                     npara(iwpara)=INT(npdz/2) + npara(iwpara-1)
                 else
                     npara(iwpara)=INT(npdz/2)
                 endif
              else
                 iwpara=iwpara-1
              endif                  
c
          elseif (paravc(1)(1:lnblnk(paravc(1))).eq.'Z_MAX()')then     

              param='-ZMAX_FILT'
              call geti_option(param,config,500,paravi,test)
              if (test.ge.1) then
                 if (iwpara.gt.1) then 
                     npara(iwpara)= test + npara(iwpara-1)
                 else
                     npara(iwpara)= test
                 endif
              else
                 iwpara=iwpara-1
              endif                  
c
          elseif (paravc(1)(1:lnblnk(paravc(1))).eq.'STRING_INPUT') then
 
              if (iwpara.gt.1) then 
                 npara(iwpara)=conti+npara(iwpara-1)
              else
                 npara(iwpara)=conti
              endif
c          
          else

              if (iwpara.gt.1) then 
                 npara(iwpara)=1+npara(iwpara-1)
              else
                 npara(iwpara)=1
              endif
ccccccccccccccccccccccc
c   check the LIB_PHYS estimators to be measured:saved in nppara
              do i = 1, 50
                if (i.le.9) write(ich,'(i1.1)') i  
                if (i.ge.10) write(ich,'(i2.2)') i  
                str1= "PHYS_PARA"// ich(1:lnblnk(ich)) //"_MED"
                str2= "PHYS_PARA"// ich(1:lnblnk(ich)) //"_INF"
                str3= "PHYS_PARA"// ich(1:lnblnk(ich)) //"_SUP"
c               
         if(paravc(1)(1:lnblnk(paravc(1))).eq.str1(1:lnblnk(str1)) .or.
     >      paravc(1)(1:lnblnk(paravc(1))).eq.str2(1:lnblnk(str2)) .or.
     >      paravc(1)(1:lnblnk(paravc(1))).eq.str3(1:lnblnk(str3)))then
                    nppara(i)=1
c                  write(UO,*) i,nppara(i),'  ',
c     >   str1(1:lnblnk(str1)),' ',str2(1:lnblnk(str2)),
c     >   '  ',str3(1:lnblnk(str3))
                endif 
              enddo
cccccccccccccccccccccc
          endif          
         endif
      enddo
 13   close(1)
c    check if output parameters defined otherwise STOP 
      if (iwpara .eq. 0) then 
        write(UO,*) 'No output parameters defined with ',
     >               file(1:lnblnk(file)), ' -->  STOP'
        STOP
      endif
      val1 = DINT(DBLE(iwpara)/5.)
      if (MOD(iwpara,5).eq.0) then
          k = IDINT(val1) 
      else
          k = IDINT(val1) + 1
      endif
      write(UO,'(A)') "#######################################"
      write(UO,'(A)') "#####        OUTPUT FORMAT        #####"
      do i = 1,k
         imin = (i-1)*5 + 1
         imax = imin + 4
         if (imax.gt.iwpara) imax = iwpara
        write(UO,'(A,1x,4(A," , "),A)') "# ",
     >  (wpara(j)(1:lnblnk(wpara(j))),j=imin,imax)
      enddo
      write(UO,'(A)') "#######################################"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
 56   write (UO,*) 'File ',file(1:lnblnk(file)),' not found -> STOP '
      stop
      end
