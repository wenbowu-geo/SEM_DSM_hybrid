!=====================================================================
!

!module constants

!  include "constants.h"

!end module constants

!=====================================================================



module struct_defined
use constants
implicit none
!******************************model for raytrace***************************
type model_1d
!Number of layers
!integer ::Nlay
!double precision,dimension(MAXLAY):: Vp,Vs,Rou,Depth
!double precision,dimension(MAXLAY):: Radius
!urn rayparameter of lth layer, Only P velocity or S velocity involved
!double precision,dimension(MAXLAY):: B,p
!integer :: Ncmb,Nicb;
end type model_1d


!***********impose traction on boundaries and predicted 1D synthetics******
type input_traction_velo
end type
!
!**************************The information needed for the integral**********
Type Variables_bound

integer                 ::ispec_bound
integer                 ::iregion  !regular(1,3,5) and unregular(2,4)
integer                 ::face_type  !5 or 6 kinds of faces(left,right,top,bottom,back,forward)
integer                 ::iele_elasitc
integer                 ::iele_acoustic
integer                 ::istart,iend,jstart,jend,kstart,kend
logical                 ::is_elastic

!real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable ::disp,traction

!double precision,dimension(NGLLX,NGLLY)   ::x,y,z
!double precision,dimension(NGLLX,NGLLY,NGLLZ)   ::sigmaxx,sigmaxy,sigmaxz,sigmayx,sigmayy,&
!                                                  sigmayz,sigmazx,sigmazy,sigmazz 
end type Variables_bound

!***************************************************************************
end module struct_defined

!***********impose traction on boundaries and predicted 1D synthetics******
module impose_1DBCS
use constants

!single variables are difend as type of double precision.
!array is type of CUSTOM_REAL.


!ipack and it_imposed corresponding to it
integer ::ipack_impose,ipack_old,it_impose

double precision::t0_this_pack

!parameters for input Green functions
double precision ::in_dt
integer  ::in_ndep_ela,in_ndep_acou
integer  ::in_maxntheta_ela,in_maxntheta_acou
integer  ::in_nsets_topsolid,in_nsets_topfluid
integer  ::in_nsets_botsolid,in_nsets_botfluid
!double precision ::in_mindep_ela,in_mindep_acou
!double precision ::in_ddep_ela,in_ddep_acou
!double precision ::in_mindist_ela,in_mindist_acou
!double precision ::in_ddist_ela,in_ddist_acou
integer  ::in_npt_eachpack

integer, dimension(:),allocatable ::ndist_ela_foridep,ndist_acou_foridep

!information for each point on each face
! id to find the corresponding traction or velocity in table for each point
integer, dimension(:,:), allocatable ::id_epsilon,id_pressure,id_velo,id_poten_dot
!integer, dimension(:,:,:), allocatable ::radiation_factor
double precision, dimension(:,:,:,:), allocatable ::zrtToxyz
real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable ::moment_zrt

!strain/pressure and velocity/potential_dot to be read
integer ::nepsilon,nvelo
integer ::npressure,npoten_dot
integer,dimension(:,:), allocatable ::impose_depPress,impose_depEpsilon
integer,dimension(:),allocatable ::idepth_used_ela,idepth_used_acou
integer :: ndepth_used_ela,ndepth_used_acou
integer,dimension(:,:,:), allocatable ::iDepDist_point
integer,dimension(:,:), allocatable ::iDepDist_elastic,iDepDist_acoustic
real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable ::impose_pressure
real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable ::impose_potential_dot
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable ::impose_epsilon
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable ::impose_velo
end module impose_1DBCS



module SEMtoTele_par
   use struct_defined
   integer ::Nele_Bound
   integer ::Nele_BoundElas
   integer ::Nele_BoundAcous
   integer ::npoints_bound
   integer ::npoints_BoundElas
   integer ::npoints_BoundAcous
   Type(Variables_bound),dimension(:),allocatable ::Bound_Info
!   real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable ::disp_bound,traction_bound
   integer, dimension(:), allocatable ::id_depth_bound
   integer, dimension(:), allocatable ::id_dist_bound

   real(kind=CUSTOM_REAL), dimension(:,:), allocatable ::normal_vect_bound
   integer ::npackage_Elas,npackage_elas_acous
   integer, dimension(:), allocatable::npoints_ipack_Elas,npoints_ipack_elas_acous
   
! for the rank=0 process
   integer ::Recv_npackage_Elas,Recv_npackage_elas_acous
   integer, dimension(:),allocatable::Recv_npoints_ipack_Elas,Recv_npoints_ipack_elas_acous

!   Type(ViaRepresent),pointer             ::RepInfo 


!used when CPML_conditions is enabled
  double precision, dimension(NDIM,NDIM) :: rotation_matrix_back_cubedsph
! global point coordinates
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore_dummy_cubedsph
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: ystore_dummy_cubedsph
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: zstore_dummy_cubedsph
!!!!!



!  USED IN MESH AND GENERATE_DATABASE
   double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,THICKNESS_BLOCK_KM, &
                   CENTER_LATITUDE_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH
   integer::NX_TOPOTAPER,NY_TOPOTAPER,NX_NOTOPO,NY_NOTOPO
   integer::TELE_IXLOW,TELE_IXHIGH,TELE_IYLOW,TELE_IYHIGH,TELE_IRTOP,TELE_IRBOT

!  USED IN SOLVER PART
   logical ::LOW_RESOLUTION
   integer ::NSTEP_BETWEEN_OUTPUTBOUND
   integer ::DECIMATE_COUPLING
   integer ::NPOINTS_PER_PACK

! USED IN MESH, GENERATE_DATABASE and SOLVER
   integer ::nxLow,nxHigh,nyLow,nyHigh,nrdown,nrtop
   integer,dimension(:),allocatable::element_xLow,element_xLowReg,element_xHigh,&
           element_xHighReg,element_yLow,element_yLowReg,element_yHigh, &
           element_yHighReg,element_rdown,element_rdownReg,element_rtop,element_rtopReg

end module SEMtoTele_Par
