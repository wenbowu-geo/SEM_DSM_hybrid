subroutine convolution()
  use coupling_SEM_DSM_par,only:nproc,myrank,npack,fft_disp_bound_new,fft_traction_bound_new,&
                          normal_vector,coord_bound,c,jacobian,NPOINTS
  use convolution_par
  use constants
  implicit none
#ifdef USE_MPI
  include "mpif.h"
#endif

  integer ::ifreq
#ifdef USE_MPI
  integer ::istat,ierr
#endif

!  call read_par_convolution()
  call allocate_array_convolution()
  call read_station()
  call calcu_rotmatrix_idGreen()

  if(myrank.eq.0) write(*,*) "rank0-Reading Green Functions and convolution starts."
    do ifreq=1,nfreq_Green
      if(myrank.eq.0.and.mod(ifreq,10).eq.0) write(*,*) "rank0-working on the ifreq=",ifreq,"(nfreq=",nfreq_Green,")"
      if(ndep_elas.gt.0) call read_GreenFunc_elastic(ifreq)
      if(ndep_acous.gt.0) call read_GreenFunc_acoustic(ifreq)
      call do_convolution(ifreq)
      if(ifreq.eq.nfreq_Green/2) write(*,*) "rank",myrank,"has done half work!"
      if(ifreq.eq.1000) write(*,*) "rank",myrank,"has done 1000 frequencies!"
    end do
  if(myrank.eq.0) write(*,*) "rank0-Convolution done."

  if(myrank.eq.0) write(*,*) "rank0-Reduce values, rotate and write seismograms."
#ifdef USE_MPI
  do istat=1,nstation
      call MPI_REDUCE(disp(1,1,istat),final_disp(1,1,istat),nfreq_green*ncomp, &
                MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(myrank.eq.0) then
         call XYZtoENZ(final_disp(1,1,istat),nfreq_Green,station_lat(istat),station_lon(istat))
         call convspc(istat,nstation,time_series_length,nfreq_green, &
               omega_imag,final_disp(1,1,istat),0,0,0,station_lat,station_lon, &
                            work_spc,work_time,station_name)
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end do
#endif

end subroutine
