!save the parameters used in the final step of DSM-SEM coupling.
subroutine save_parameter_coupling(LOCAL_PATH)
use constants
use SEMtoTele_par
use generate_databases_par, only: DT,NSTEP
use mpi
implicit none

character(len=256) LOCAL_PATH

integer ::nstep_after_decimating,nsection_each_package
character(len=256) ::final_LOCAL_PATH,clean_LOCAL_PATH


!  if( USE_OUTPUT_FILES_PATH ) then
!      final_LOCAL_PATH = OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)) //'/'
!  else
      ! suppress white spaces if any
      clean_LOCAL_PATH = adjustl(LOCAL_PATH)
      ! create full final local path
      final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) //'/'
!  endif


! open file

!open(unit=Imain_SEM_Par_Coupling,file=MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES)) &
!       //'SEM_Par_Coupling',status='unknown',action='write',form='formatted')

open(unit=Imain_SEM_Par_Coupling,file= trim(trim(final_LOCAL_PATH)//"SEM_Par_Coupling"),&
       status='unknown',action='write',form='formatted')


write(Imain_SEM_Par_Coupling,*)  'npackage_SEM = ',npackages_total_coupling
write(Imain_SEM_Par_Coupling,*)  'npoints_per_pack_SEM = ', NPOINTS_PER_PACK
write(Imain_SEM_Par_Coupling,*)  'nstep_each_section_SEM = ', NSTEP_BETWEEN_OUTPUTBOUND
nstep_after_decimating=int(NSTEP/DECIMATE_COUPLING)
write(Imain_SEM_Par_Coupling,*)  'total_nstep_SEM = ', nstep_after_decimating
!Take the nearest integer on the right hand side
nsection_each_package=nint((nstep_after_decimating+TINYVAL)/NSTEP_BETWEEN_OUTPUTBOUND)
write(Imain_SEM_Par_Coupling,*)  'nsection_SEM = ', nsection_each_package
write(Imain_SEM_Par_Coupling,*)  'deltat_SEM = ', DT*DECIMATE_COUPLING
close(Imain_SEM_Par_Coupling)

end subroutine save_parameter_coupling

subroutine save_package_list(myrank,LOCAL_PATH,nproc)
  use SEMtoTele_par
  use mpi
  use constants, only:DEBUG_COUPLING
  implicit none
!  include 'mpif.h'

  integer ::myrank,nproc
  character(len=256) LOCAL_PATH

  integer ::iproc,ipack_iproc0,ipack,ipack_local,ier
  integer, dimension(MPI_STATUS_SIZE)  :: request_mpi_status
  character(len=256) ::package_list_file,final_LOCAL_PATH,clean_LOCAL_PATH


!  if( USE_OUTPUT_FILES_PATH ) then
!      final_LOCAL_PATH = OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)) //'/'
!  else
      ! suppress white spaces if any
      clean_LOCAL_PATH = adjustl(LOCAL_PATH)
      ! create full final local path
      final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) //'/'
!  endif

  package_list_file=trim(final_LOCAL_PATH)//"package_list"
  if(myrank.eq.0) then
    open(unit=121,file=trim(package_list_file),status='unknown',&
         action='write',form='formatted',iostat=ier)
    
    ipack=0
    do ipack_iproc0=1,npackage_elas_acous
      ipack=ipack+1
      write(121,*)ipack,0,ipack_iproc0,npoints_ipack_elas_acous(ipack_iproc0)
    end do
    do iproc=1,nproc-1
       call MPI_RECV(Recv_npackage_elas_acous,1,MPI_INTEGER,iproc,0,MPI_COMM_WORLD, &
                      request_mpi_status,ier)
       if(DEBUG_COUPLING) print *,'Recv_npackage_elas_acous', Recv_npackage_elas_acous
       allocate(Recv_npoints_ipack_elas_acous(Recv_npackage_elas_acous))
       call MPI_RECV(Recv_npoints_ipack_elas_acous,Recv_npackage_elas_acous,MPI_INTEGER, &
                      iproc,1,MPI_COMM_WORLD,request_mpi_status,ier)
       do ipack_local=1,Recv_npackage_elas_acous
          ipack=ipack+1
          write(121,*)ipack,iproc,ipack_local,Recv_npoints_ipack_elas_acous(ipack_local)
       end do
       npackages_total_coupling=ipack
       deallocate(Recv_npoints_ipack_elas_acous)
    end do

  else
    call MPI_SEND(npackage_elas_acous,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ier)
    call MPI_SEND(npoints_ipack_elas_acous,npackage_elas_acous,MPI_INTEGER,0,1, &
                      MPI_COMM_WORLD,ier)
  end if
  if(myrank.eq.0) close(121)
  

end subroutine save_package_list

