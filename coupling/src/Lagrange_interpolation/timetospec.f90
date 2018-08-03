program coupling_SEM_DSM
    use coupling_SEM_DSM_par

    implicit none
#ifdef USE_MPI
    include "mpif.h"
#endif

    integer ::ipackage,ipoint
#ifdef USE_MPI
    integer ::ier
#endif


    call initialize()
    call read_para()
    call partition_job()
    call allocate_array1()
    call read_package_id()
    call allocate_array2()
    call read_par_convolution()

    if(myrank.eq.0) write(*,*) "Start reading boundary displacement and traction."
    do ipackage=1,npack
!      print *,'read ipackage',ipackage
      call read_ipack(ipack_start+ipackage-1,ipackage)
      do ipoint=1,npoint_pack(ipack_start+ipackage-1)
        call resample(ipoint,ipackage)
        call fft(ipoint-1+ipoint_start(ipackage))
      end do
    end do
    if(myrank.eq.0) write(*,*) "Reading boundary displacement and traction done."
    call get_boundinfo()
    call deallocate_array()
!    call MPI_BARRIER(MPI_COMM_WORLD,ier)
    if(myrank.eq.0) write(*,*) "Convolution starts."
    call convolution()
    if(myrank.eq.0) write(*,*) "Done!"
#ifdef USE_MPI
   call MPI_FINALIZE(ier)
#endif


end program
