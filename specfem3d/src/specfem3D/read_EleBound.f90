subroutine read_EleBound()
   use specfem_par
   use specfem_par_acoustic
   use specfem_par_elastic
   use specfem_par_poroelastic
   use SEMtoTele_par

   implicit none
!   include "mpif.h"
!   character(len=256) proname,file_name
   integer ::i,ier

   call create_name_database(prname,myrank,LOCAL_PATH)
   open(unit=175,file=prname(1:len_trim(prname))//'TeleEle_xLow.info',status='old',&
   action='read',form='formatted',iostat=ier)
   read(175,*) nxLow
   if(nxLow.gt.0) then
        allocate(element_xLow(nxLow))
        allocate(element_xLowReg(nxLow))
   end if
   do i=1,nxLow
        read(175,*)element_xLow(i),element_xLowReg(i)
   end do
   close(175)

   open(unit=175,file=prname(1:len_trim(prname))//'TeleEle_xHigh.info',status='old',&
      action='read',form='formatted',iostat=ier)
   read(175,*) nxHigh
   if(nxHigh.gt.0) then
        allocate(element_xHigh(nxHigh))
        allocate(element_xHighReg(nxHigh))
   end if
   do i=1,nxHigh
        read(175,*)element_xHigh(i),element_xHighReg(i)
   end do
   close(175)


   open(unit=175,file=prname(1:len_trim(prname))//'TeleEle_yLow.info',status='old',&
      action='read',form='formatted',iostat=ier)
  read(175,*) nyLow
  if(nyLow.gt.0) then
     allocate(element_yLow(nyLow))
     allocate(element_yLowReg(nyLow))
  end if
  do i=1,nyLow
        read(175,*)element_yLow(i),element_yLowReg(i)
  end do
  close(175)


  open(unit=175,file=prname(1:len_trim(prname))//'TeleEle_yHigh.info',status='old',&
     action='read',form='formatted',iostat=ier)
   read(175,*) nyHigh
   if(nyHigh.gt.0) then
         allocate(element_yHigh(nyHigh))
         allocate(element_yHighReg(nyHigh))
   end if
   do i=1,nyHigh
         read(175,*)element_yHigh(i),element_yHighReg(i)
   end do
   close(175)
   
  open(unit=175,file=prname(1:len_trim(prname))//'TeleEle_rTop.info',status='old',&
     action='read',form='formatted',iostat=ier)
  read(175,*) nrtop
  if(nrtop.gt.0) then
        allocate(element_rtop(nrtop))
        allocate(element_rtopReg(nrtop))
  end if
  do i=1,nrtop
        read(175,*)element_rtop(i),element_rtopReg(i)
  end do
  close(175)

  open(unit=175,file=prname(1:len_trim(prname))//'TeleEle_rBot.info',status='old',&
     action='read',form='formatted',iostat=ier)
  read(175,*) nrdown
  if(nrdown.gt.0) then
        allocate(element_rdown(nrdown))
        allocate(element_rdownReg(nrdown))
  end if
  do i=1,nrdown
        read(175,*)element_rdown(i),element_rdownReg(i)
  end do
  close(175)
end subroutine read_EleBound
