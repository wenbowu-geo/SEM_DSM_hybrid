subroutine save_SEMtoTeleEle(TeleEle_nxLow,TeleEle_xLow,TeleEle_xLowReg,TeleEle_nxHigh,TeleEle_xHigh, &
                TeleEle_xHighReg,TeleEle_nyLow,TeleEle_yLow,TeleEle_yLowReg,TeleEle_nyHigh,TeleEle_yHigh,&
                TeleEle_yHighReg,TeleEle_nrTop,TeleEle_rTop,TeleEle_rTopReg,&
                TeleEle_nrBot,TeleEle_rBot,TeleEle_rBotReg,LOCAL_PATH,iproc)

  use constants
  implicit none
!  include "constants.h"
  ! name of the database files
  integer ::iproc
  character(len=256) proname,file_name,LOCAL_PATH,clean_LOCAL_PATH
  integer ::TeleEle_nxLow,TeleEle_nxHigh,TeleEle_nyLow,TeleEle_nyHigh,TeleEle_nrTop,TeleEle_nrBot
  integer,dimension(TeleEle_nxLow)   ::TeleEle_xLow,TeleEle_xLowReg
  integer,dimension(TeleEle_nxHigh)  ::TeleEle_xHigh,TeleEle_xHighReg
  integer,dimension(TeleEle_nyLow)   ::TeleEle_yLow,TeleEle_yLowReg
  integer,dimension(TeleEle_nyHigh)  ::TeleEle_yHigh,TeleEle_yHighReg
  integer,dimension(TeleEle_nrTop)      ::TeleEle_rTop,TeleEle_rTopReg
  integer,dimension(TeleEle_nrBot)      ::TeleEle_rBot,TeleEle_rBotReg

  integer i


  write(proname,"('/proc',i6.6,'_')") iproc
  clean_LOCAL_PATH = adjustl(LOCAL_PATH)
  file_name = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH))//proname(1:len_trim(proname))//'TeleEle_xLow.info'
!  print *,'file',iproc,proname(1:len_trim(proname))//'TeleEle_xLow.info'

  open(unit=16,file=file_name(1:len_trim(file_name)),status='unknown',action='write',form='formatted')
  write(16,*) TeleEle_nxLow
  do i=1,TeleEle_nxLow
     write(16,*)TeleEle_xLow(i),TeleEle_xLowReg(i)
  end do
  close(16)
  
  file_name =  clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH))//proname(1:len_trim(proname))//'TeleEle_xHigh.info'

  open(unit=16,file=file_name(1:len_trim(file_name)),status='unknown',action='write',form='formatted')
  write(16,*) TeleEle_nxHigh
  do i=1,TeleEle_nxHigh
     write(16,*)TeleEle_xHigh(i),TeleEle_xHighReg(i)
  end do
  close(16)

  file_name = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH))//proname(1:len_trim(proname))//'TeleEle_yLow.info'
  open(unit=16,file=file_name(1:len_trim(file_name)),status='unknown',action='write',form='formatted')
  write(16,*) TeleEle_nyLow
  do i=1,TeleEle_nyLow
     write(16,*)TeleEle_yLow(i),TeleEle_yLowReg(i)
  end do
  close(16)

   file_name = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH))//proname(1:len_trim(proname))//'TeleEle_yHigh.info'
  open(unit=16,file=file_name(1:len_trim(file_name)),status='unknown',action='write',form='formatted')
  write(16,*) TeleEle_nyHigh
  do i=1,TeleEle_nyHigh
     write(16,*)TeleEle_yHigh(i),TeleEle_yHighReg(i)
  end do
  close(16)

  file_name =  clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH))//proname(1:len_trim(proname))//'TeleEle_rTop.info'
  open(unit=16,file=file_name(1:len_trim(file_name)),status='unknown',action='write',form='formatted')
  if(dabs(R_TOP_BOUND-R_EARTH_SURF).lt.TINYVAL) then
!to be fixed. It should be correct, but double check.
! The top boundary is free surface, which has no contrinutions in coupling.
    write(16,*) 0
  else
    write(16,*) TeleEle_nrTop
    do i=1,TeleEle_nrTop
       write(16,*)TeleEle_rTop(i),TeleEle_rTopReg(i)
    end do
    close(16)
  end if

  file_name =  clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH))//proname(1:len_trim(proname))//'TeleEle_rBot.info'
  open(unit=16,file=file_name(1:len_trim(file_name)),status='unknown',action='write',form='formatted')
  write(16,*) TeleEle_nrBot
  do i=1,TeleEle_nrBot
     write(16,*)TeleEle_rBot(i),TeleEle_rBotReg(i)
  end do
  close(16)



end subroutine

