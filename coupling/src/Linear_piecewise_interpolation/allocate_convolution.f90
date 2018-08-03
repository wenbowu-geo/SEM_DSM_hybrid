subroutine allocate_array_convolution()
use convolution_par
use coupling_SEM_DSM_par,only:npoints,myrank
use constants

implicit none

allocate(disp(ncomp,nfreq_Green,nstation))
!if(myrank.eq.0) then
  allocate(final_disp(ncomp,nfreq_Green,nstation))
!end if
allocate(station_lat(nstation))
allocate(station_lon(nstation))
allocate(station_name(nstation))
allocate(station_coord(ncomp,nstation))


allocate(disp_Green_elas(ncomp,nstat_Green_elas,ncomp))
allocate(epsilon_Green(nstat_Green_elas,ncomp_stress,ncomp))
allocate(disp_Green_acous(ncomp,nstat_Green_acous,ncomp))
allocate(pressure_Green(nstat_Green_acous,ncomp))

allocate(rot_matrix_bound(ncomp,ncomp,npoints,nstation))
allocate(rot_matrix_station(ncomp,ncomp,npoints,nstation))
!allocate(gcarc(npoints,nstation))
allocate(id_Green(npoints,nstation))
allocate(igcarc_remained(npoints,nstation))

allocate(work_time(32*nfreq_Green))
allocate(work_spc(16*nfreq_Green))

disp(:,:,:)=cmplx(0.d0)
if(myrank.eq.0) then
   final_disp(:,:,:)=cmplx(0.d0)
end if
disp_Green_elas(:,:,:)=cmplx(0.d0)
epsilon_Green(:,:,:)=cmplx(0.d0)
disp_Green_elas(:,:,:)=cmplx(0.d0)
pressure_Green(:,:)=cmplx(0.d0)
rot_matrix_bound(:,:,:,:)=0.d0
rot_matrix_station(:,:,:,:)=0.d0
!gcarc(:,:)=0.d0
id_Green(:,:)=0
station_coord(:,:)=0.0


work_time(:)=0.d0
work_spc(:)=cmplx(0.d0)
end subroutine allocate_array_convolution
