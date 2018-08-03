subroutine allocate_array1()
use coupling_SEM_DSM_par

implicit none

allocate(global_pack_id(npackage_SEM))
allocate(local_pack_id(npackage_SEM))
allocate(in_iproc(npackage_SEM))
allocate(NPOINT_PACK(npackage_SEM))
allocate(ipoint_start(npack))

global_pack_id(:)=0
local_pack_id(:)=0
in_iproc(:)=0
NPOINT_PACK(:)=0
ipoint_start(:)=0

end subroutine allocate_array1

subroutine allocate_array2()
use coupling_SEM_DSM_par
use constants

implicit none

integer ::array_size
array_size=max(total_nstep_SEM*ncomp*NPOINTS,double_nfreq_DSM*ncomp*NPOINTS)
!print *,"max_arraysize_allowd is",max_array_size,'and the measured size is ',array_size,npackage_SEM
if(array_size.gt.max_array_size) then
   call exit_MPI('the size of array is too large')
end if

if(myrank.eq.0) print *,"allocate_array2 starts ..."
allocate(disp_bound(total_nstep_SEM,ncomp,npoints_per_pack_SEM))
allocate(traction_bound(total_nstep_SEM,ncomp,npoints_per_pack_SEM))
!not used
!allocate(disp_bound_new(double_nfreq_DSM,ncomp,npoints_per_pack_SEM))
!allocate(traction_bound_new(double_nfreq_DSM,ncomp,npoints_per_pack_SEM))
allocate(disp_bound_fine(npts_fine,ncomp))
allocate(traction_bound_fine(npts_fine,ncomp))


allocate(fft_disp_bound_new(double_nfreq_DSM,ncomp,NPOINTS))
allocate(fft_traction_bound_new(double_nfreq_DSM,ncomp,NPOINTS))
allocate(fft_disp_bound_fine(npts_fine,ncomp))
allocate(fft_traction_bound_fine(npts_fine,ncomp))

allocate(id_depth(NPOINTS))
allocate(is_elastic(NPOINTS))
allocate(is_acoustic(NPOINTS))
allocate(normal_vector(ncomp,NPOINTS))
allocate(jacobian(NPOINTS))
allocate(coord_bound(ncomp,NPOINTS))
allocate(c(N_ELAS_COEF,NPOINTS))

!allocate(new_local_global(npoint))

disp_bound(:,:,:)=0.d0
traction_bound(:,:,:)=0.d0
!disp_bound_new(:,:,:)=0.d0
!traction_bound_new(:,:,:)=0.d0
disp_bound_fine(:,:)=0.d0
traction_bound_fine(:,:)=0.d0
fft_disp_bound_new(:,:,:)=0.d0
fft_traction_bound_new(:,:,:)=0.d0
fft_disp_bound_fine(:,:)=0.d0
fft_traction_bound_fine(:,:)=0.d0

!new_local_global(:)=0
normal_vector(:,:)=0.d0
jacobian(:)=0.d0
coord_bound(:,:)=0.d0
c(:,:)=0.d0
if(myrank.eq.0) print *,"allocate_array2 finishs!"

end subroutine allocate_array2

subroutine deallocate_array()

use coupling_SEM_DSM_par
deallocate(disp_bound)
deallocate(traction_bound)
!deallocate(disp_bound_new)
!deallocate(traction_bound_new)
deallocate(fft_disp_bound_fine)
deallocate(fft_traction_bound_fine)
deallocate(disp_bound_fine)
deallocate(traction_bound_fine)
deallocate(global_pack_id)
deallocate(local_pack_id)
deallocate(in_iproc)
deallocate(NPOINT_PACK)
deallocate(ipoint_start)

end subroutine

