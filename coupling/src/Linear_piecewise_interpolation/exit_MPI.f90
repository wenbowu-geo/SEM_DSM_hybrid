subroutine exit_MPI(error_msg)

  implicit none
#ifdef USE_MPI
  ! standard include of the MPI library
  include "mpif.h"
#endif

  ! identifier for error message file
  integer, parameter :: IERROR = 30

  character(len=*) error_msg

  integer ier

  ier = 0

  ! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI... proc '

  ! stop all the MPI processes, and exit
#ifdef USE_MPI
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
#endif

  stop 'error, program ended in exit_MPI'

end subroutine exit_MPI

