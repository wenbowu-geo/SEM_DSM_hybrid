subroutine read_para()
   use coupling_SEM_DSM_par
   use convolution_par
   use constants
   implicit none

! other variables
   double precision ::npow,ratio_deltat
   integer ::n_refine
   character(len=80) ::para_file,para_file_SEM_Par,GREEN_par_file
   integer,parameter ::IIN_SEM_Par_Coupling=14
   integer,parameter ::IIN_Green_Par=15
   integer ::ios
   integer ::ier

   para_file="DATA/Par_file"
   open(unit=IIN,file=trim(para_file),action='read',form="formatted",status="old",iostat=ios)
   if(ios /= 0) then
      print *,'error:could not open parameter file'
      call exit_mpi('error opening parameter file')
   end if

   call read_value_integer(IIN,IGNORE_JUNK,npoints_taper_SEM,'coupling.npoints_taper_SEM',ier)
   if(ier /= 0) stop 'Error reading coupling parameter npoints_taper_SEM'

   call read_value_logical(IIN,IGNORE_JUNK,do_coupling_RTZ(1),'coupling.vertical_component',ier)
   if(ier /= 0) stop 'Error reading coupling parameter vertical_component'

   call read_value_logical(IIN,IGNORE_JUNK,do_coupling_RTZ(2),'coupling.radial_component',ier)
   if(ier /= 0) stop 'Error reading coupling parameter radial_component'

   call read_value_logical(IIN,IGNORE_JUNK,do_coupling_RTZ(3),'coupling.transverse_component',ier)
   if(ier /= 0) stop 'Error reading coupling parameter transverse_component'

   close (IIN)



   para_file_SEM_Par="SEM_input/SEM_Par_Coupling"
   open(unit=IIN_SEM_Par_Coupling,file=trim(para_file_SEM_Par),action='read',&
        form="formatted",status="old",iostat=ios)
   if(ios /= 0) then
      print *,'error:could not open SEM_input/SEM_Par_Coupling'
      call exit_mpi('error opening SEM_input/SEM_Par_Coupling')
   end if

   call read_value_integer(IIN_SEM_Par_Coupling,IGNORE_JUNK,npackage_SEM,'SEM_Par.npackage_SEM',ier)
   if(ier /= 0) stop 'Error reading SEM parameter npackage_SEM'

   call read_value_integer(IIN_SEM_Par_Coupling,IGNORE_JUNK,npoints_per_pack_SEM,'SEM_Par.npoints_per_pack_SEM',ier)
   if(ier /= 0) stop 'Error reading SEM parameter npoints_per_pack_SEM'

   call read_value_integer(IIN_SEM_Par_Coupling,IGNORE_JUNK,nstep_each_section_SEM,'SEM_Par.nstep_each_section_SEM',ier)
   if(ier /= 0) stop 'Error reading SEM parameter nstep_each_section_SEM'

   call read_value_integer(IIN_SEM_Par_Coupling,IGNORE_JUNK,total_nstep_SEM,'SEM_Par.total_nstep_SEM',ier)
   if(ier /= 0) stop 'Error reading SEM parameter total_nstep_SEM'

   call read_value_integer(IIN_SEM_Par_Coupling,IGNORE_JUNK,nsection_SEM,'SEM_Par.nsection_SEM',ier)
   if(ier /= 0) stop 'Error reading SEM parameter nsection_SEM'

   call read_value_dble_precision(IIN_SEM_Par_Coupling,IGNORE_JUNK,deltat_SEM,'SEM_Par.deltat_SEM',ier)
   if(ier /= 0) stop 'Error reading SEM parameter deltat_SEM'
   
   close(IIN_SEM_Par_Coupling)


   Green_par_file="./DSM_input/Green_Par"
   open(unit=IIN_Green_Par,file=trim(Green_par_file),status='unknown',form='formatted',iostat=ios)
   if(ios /= 0) stop 'Error in openning file Green_Par'
   read(IIN_Green_Par,*)nfreq_Green
   read(IIN_Green_Par,*)min_dep,ddep,ndep_acous
   read(IIN_Green_Par,*)min_dep,ddep,ndep_elas
   read(IIN_Green_Par,*)min_theta,dtheta,ntheta
   read(IIN_Green_Par,*)time_series_length
   read(IIN_Green_Par,*)omega_imag
   close(IIN_Green_Par)
   nstat_Green_elas=ndep_elas*ntheta
   nstat_Green_acous=ndep_acous*ntheta
   min_dep=min_dep*km
   ddep=ddep*km




   double_nfreq_DSM=2*nfreq_Green
   deltat_DSM = time_series_length/double_nfreq_DSM
   npow=double_nfreq_DSM/2.0
   do while(npow.gt.1.d0)
      npow=npow/2.d0
   end do
   if(dabs(npow-1.d0)>TINY) then 
     call exit_MPI('error of double_nfreq_DSM, double_nfreq_DSM should be power of 2')
   end if

   ratio_deltat=deltat_DSM/deltat_SEM
   n_refine=1
   if(ratio_deltat.lt.1.d0+TINY) then
     call exit_MPI('error of deltat, sampling for the boundary disp is too large')
   else
     do while(ratio_deltat.gt.1.d0)
        ratio_deltat=ratio_deltat/2.d0
        n_refine=n_refine*2
     end do
   end if

!to be fixed
!   n_refine=n_refine/4
   n_refine=n_refine/2

!   print *,'n_refine',n_refine,deltat_refine,deltat_DSM
   npts_fine=double_nfreq_DSM*n_refine
   deltat_refine=deltat_DSM/n_refine
   !print *,'n_refine',n_refine,deltat_refine,deltat_DSM

   !if(n_refine.le.2) then
   !    deltat_refine=deltat_DSM
   !    npts_fine=double_nfreq_DSM
   !else
   !    deltat_refine=deltat_DSM/n_refine*2
   !    npts_fine=double_nfreq_DSM*n_refine/2
   !end if
end subroutine read_para



subroutine read_package_id()
    use coupling_SEM_DSM_par
!other variables
    integer ::ipack,ios
    character(len=80) ::package_file

   package_file="SEM_input/package_list"
   open(unit=10,file=trim(package_file),action='read',form="formatted",status="old",iostat=ios)
   if(ios /= 0) then
      print *,'error:could not open station_list file'
      call exit_mpi('error opening station_list file')
   end if
   
   do ipack=1,npackage_SEM
       read(10,*) global_pack_id(ipack),in_iproc(ipack),local_pack_id(ipack),npoint_pack(ipack)
!       read(*) global_pack_id(ipack),in_iproc(ipack),local_pack_id(ipack),npoint_pack(ipack)
   end do
   close(10)


   npoints=0
   do ipack=1,npack
     if(ipack.eq.1) then
      ipoint_start(ipack)=1
     else 
      ipoint_start(ipack)=ipoint_start(ipack-1)+npoint_pack(ipack-1+ipack_start-1)
     end if
   end do
   npoints=ipoint_start(npack)+npoint_pack(npack+ipack_start-1)-1

end subroutine read_package_id


subroutine read_ipack(ipack,ipack_local)
    use coupling_SEM_DSM_par
    implicit none
#ifdef USE_MPI
  ! standard include of the MPI library
  include "mpif.h"
#endif

    integer ::ipack,ipack_local

!other variables
   integer ::iblock_time,ipoint,nstep_this_block,it_this_block,it
!integer   ::icomp
   character(len=80) ::disp_file,traction_file
   integer ::ios
   integer ::ier

!   double precision ::nstation_temp



!**********************************************************************
!**************************read the displacement file******************
   write(disp_file,"('./SEM_input/iproc',i4.4,'/disp_pack',i6.6)")in_iproc(ipack),local_pack_id(ipack)
   open(unit=30,file=trim(disp_file),action='read',access='stream',form='unformatted',status="old",iostat=ios)
!   open(unit=30,file=trim(disp_file),action='read',form="formatted",status="old",iostat=ios)
   if(ios /= 0) then
        print *,'displace file',disp_file
        call exit_mpi('error opening disp file')
   end if

!   read(30) nstation_temp
!   if(nstation_temp.ne.npoint(ipack)) then
!        call exit_mpi('error reading nstation')
!   end if
!   call MPI_BARRIER(MPI_COMM_WORLD,ier)
   do iblock_time=1,nsection_SEM
    do ipoint=1,npoint_pack(ipack)
     if(iblock_time.lt.nsection_SEM) then
        nstep_this_block=nstep_each_section_SEM
     else if(iblock_time.eq.nsection_SEM) then
        nstep_this_block=total_nstep_SEM-nstep_each_section_SEM*(nsection_SEM-1)
     end if
     do it_this_block=1,nstep_this_block
        if(ipoint.gt.npoints_per_pack_SEM)then
               stop 'Error reading'
        end if
        it=it_this_block+nstep_each_section_SEM*(iblock_time-1)
        read(30) disp_bound(it,:,ipoint)
     end do
    end do
   end do
   close(30)
!   call MPI_BARRIER(MPI_COMM_WORLD,ier)

!**********************************************************************
!**************************read the stress file******************
   write(traction_file,"('./SEM_input/iproc',i4.4,'/traction_pack',i6.6)")in_iproc(ipack),local_pack_id(ipack)
   open(unit=40,file=trim(traction_file),action='read',access='stream',form='unformatted',status="old",iostat=ios)
!   open(unit=40,file=trim(traction_file),action='read',form="formatted",status="old",iostat=ios)
   if(ios /= 0) then
           call exit_mpi('error opening traction file')
   end if

!   read(40) nstation_temp
   do iblock_time=1,nsection_SEM
    do ipoint=1,npoint_pack(ipack)
     if(iblock_time.lt.nsection_SEM) then
        nstep_this_block=nstep_each_section_SEM
     else if(iblock_time.eq.nsection_SEM) then
        nstep_this_block=total_nstep_SEM-nstep_each_section_SEM*(nsection_SEM-1)
     end if
     do it_this_block=1,nstep_this_block
       if(ipoint.gt.npoints_per_pack_SEM) stop 'Error reading'
!       read(40) traction_bound(it,:,ipoint-1+ipoint_start(ipack_local))
        it=it_this_block+nstep_each_section_SEM*(iblock_time-1)
!        read(40,*) traction_bound(it,:,ipoint)
        read(40) traction_bound(it,:,ipoint)
     end do
    end do
   end do
   close(40)


end subroutine read_ipack



subroutine read_number_stations()
   use convolution_par
   use constants, only:km
   implicit none

!other variables
   character(len=80) ::station_file
   integer ::ios

  station_file="./DATA/STATION"
  open(unit=40,file=trim(station_file),action='read',form="formatted",status="unknown",iostat=ios)
  if(ios /= 0) then
    call exit_MPI('Error in openning file DATA/STATION')
  end if
  read(40,*)nstation
  close(40)


end subroutine read_number_stations

subroutine read_station()
   use convolution_par
   implicit none
!other variables
  integer ::ios,nstation_temp,ista
  character(len=80) ::station_file

  station_file="./DATA/STATION"
  open(unit=40,file=trim(station_file),action='read',form="formatted",status="unknown",iostat=ios)
  if(ios /= 0) then
    call exit_MPI('Error in openning file DATA/STATION')
  end if
  read(40,*) nstation_temp
  do ista=1,nstation
     read(40,*) station_name(ista),station_lon(ista),station_lat(ista)
  end do
  close(40)
end subroutine read_station

