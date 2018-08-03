subroutine save_package_list(myrank,LOCAL_PATH,nproc)
  use SEMtoTele_par
  use mpi
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
       print *,'Recv_npackage_elas_acous', Recv_npackage_elas_acous
       allocate(Recv_npoints_ipack_elas_acous(Recv_npackage_elas_acous))
       call MPI_RECV(Recv_npoints_ipack_elas_acous,Recv_npackage_elas_acous,MPI_INTEGER, &
                      iproc,1,MPI_COMM_WORLD,request_mpi_status,ier)
       do ipack_local=1,Recv_npackage_elas_acous
          ipack=ipack+1
          write(121,*)ipack,iproc,ipack_local,Recv_npoints_ipack_elas_acous(ipack_local)
       end do
       deallocate(Recv_npoints_ipack_elas_acous)
    end do

  else
    call MPI_SEND(npackage_elas_acous,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ier)
    call MPI_SEND(npoints_ipack_elas_acous,npackage_elas_acous,MPI_INTEGER,0,1, &
                      MPI_COMM_WORLD,ier)
  end if
  if(myrank.eq.0) close(121)
  

end subroutine save_package_list
