module constants
   include "constants.h"
end module constants

module coupling_SEM_DSM_par
  integer ::nproc,myrank
  
  integer ::npackage_SEM,npack,npoints_per_pack_SEM
  integer ::ipack_start,ipack_end
  integer ::total_nstep_SEM,double_nfreq_DSM,npoints_taper_SEM,npts_fine
  integer ::nsection_SEM,nstep_each_section_SEM
  double precision ::deltat_SEM,deltat_DSM,deltat_refine
  logical,dimension(3) ::do_coupling_RTZ
 
  
  integer, dimension(:),allocatable ::in_iproc,global_pack_id,local_pack_id
  integer, dimension(:),allocatable ::npoint_pack
  integer                           ::npoints
  integer, dimension(:),allocatable  ::ipoint_start
  double precision,dimension(:,:,:), allocatable ::disp_bound,traction_bound
  double precision,dimension(:,:,:), allocatable ::disp_bound_new,traction_bound_new
 
! temperary array for the resampling and fft calculation for each point. we want get get an good sampling rate 
!which makes the number of interpolation points not very high(eg. 8-16 points).
! This rate should also be npower(2)*npts_new. Then we can directly copy this value to fft_disp_bound_new.
  double precision,dimension(:,:),allocatable ::disp_bound_fine,traction_bound_fine


  !integer ::nfreq
  complex(kind=8),dimension(:,:,:),allocatable ::fft_disp_bound_new,fft_traction_bound_new
  complex(kind=8),dimension(:,:),allocatable::fft_disp_bound_fine,fft_traction_bound_fine

  integer, dimension(:), allocatable ::id_depth
  logical, dimension(:), allocatable ::is_elastic,is_acoustic
  double precision,dimension(:,:), allocatable ::normal_vector,coord_bound,c
  double precision,dimension(:), allocatable   ::jacobian
  
 
end module coupling_SEM_DSM_par
