!=====================================================================
!WENBO
!This subroutine is similar with original code's subroutine 'get_flags_boundaries',
!Just modifying it for the cubed-sphere projection.
!=====================================================================

  subroutine get_flags_boundaries_cubed(nspec,iproc_xi,iproc_eta,ispec,idoubling, &
             xelm,yelm,zelm,iboun,iMPIcut_xi,iMPIcut_eta, &
             NPROC_XI,NPROC_ETA,ANGULAR_WIDTH_XI_IN_DEGREES, &
             ANGULAR_WIDTH_ETA_IN_DEGREES,Z_DEPTH_BLOCK,IFLAG_ONE_LAYER_TOPOGRAPHY)

  use constants
  implicit none

!  include "constants.h"

  integer nspec
  integer ispec,idoubling
  integer NPROC_XI,NPROC_ETA

  double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,Z_DEPTH_BLOCK
  integer IFLAG_ONE_LAYER_TOPOGRAPHY

  logical iboun(6,nspec)
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)


! use iproc_xi and iproc_eta to determine MPI cut planes along xi and eta
  integer iproc_xi,iproc_eta

  double precision target,sizeslice,TOLERANCE_DEGREES,TOLERANCE_METERS
  double precision xelm(8),yelm(8),zelm(8)

! compute geometrical tolerance small compared to size of model to detect edges
  TOLERANCE_DEGREES = dabs(ANGULAR_WIDTH_XI_IN_DEGREES) / 100000.
  TOLERANCE_METERS = dabs(Z_DEPTH_BLOCK)/100000
! ****************************************************
!     determine if the element falls on a boundary
! ****************************************************

  iboun(:,ispec)=.false.

! on boundary 1: x=xmin
  target= TOLERANCE_DEGREES
  if(xelm(1)<target .and. xelm(4)<target .and. xelm(5)<target .and. xelm(8)<target) iboun(1,ispec)=.true.

! on boundary 2: xmax
  target= ANGULAR_WIDTH_XI_IN_DEGREES - TOLERANCE_DEGREES
  if(xelm(2)>target .and. xelm(3)>target .and. xelm(6)>target .and. xelm(7)>target) iboun(2,ispec)=.true.

! on boundary 3: ymin
  target= TOLERANCE_DEGREES
  if(yelm(1)<target .and. yelm(2)<target .and. yelm(5)<target .and. yelm(6)<target) iboun(3,ispec)=.true.

! on boundary 4: ymax
  target= ANGULAR_WIDTH_ETA_IN_DEGREES - TOLERANCE_DEGREES
  if(yelm(3)>target .and. yelm(4)>target .and. yelm(7)>target .and. yelm(8)>target) iboun(4,ispec)=.true.

! on boundary 5: bottom
  target = Z_DEPTH_BLOCK + TOLERANCE_METERS
  if(zelm(1)<target .and. zelm(2)<target .and. zelm(3)<target .and. zelm(4)<target) iboun(5,ispec)=.true.
!WENBO
! on boundary 6: top
  if(idoubling == IFLAG_ONE_LAYER_TOPOGRAPHY) iboun(6,ispec)=.true.
!WENBO


! *******************************************************************
!     determine if the element falls on an MPI cut plane along xi
! *******************************************************************

! detect the MPI cut planes along xi in the cubed sphere

  iMPIcut_xi(:,ispec)=.false.

! angular size of a slice along xi
  sizeslice = (ANGULAR_WIDTH_XI_IN_DEGREES) / NPROC_XI

! left cut-plane in the current slice along X = constant (Xmin of this slice)
! and add geometrical tolerance

  target = iproc_xi*sizeslice + TOLERANCE_DEGREES
  if(xelm(1)<target .and. xelm(4)<target .and. xelm(5)<target .and. xelm(8)<target) &
    iMPIcut_xi(1,ispec)=.true.

! right cut-plane in the current slice along X = constant (Xmax of this slice)
! and add geometrical tolerance

  target =(iproc_xi+1)*sizeslice - TOLERANCE_DEGREES
  if(xelm(2)>target .and. xelm(3)>target .and. xelm(6)>target .and. xelm(7)>target) &
    iMPIcut_xi(2,ispec)=.true.

! ********************************************************************
!     determine if the element falls on an MPI cut plane along eta
! ********************************************************************

  iMPIcut_eta(:,ispec)=.false.

! angular size of a slice along eta
  sizeslice = (ANGULAR_WIDTH_ETA_IN_DEGREES) / NPROC_ETA

! left cut-plane in the current slice along Y = constant (Ymin of this slice)
! and add geometrical tolerance

  target = iproc_eta*sizeslice + TOLERANCE_DEGREES
  if(yelm(1)<target .and. yelm(2)<target .and. yelm(5)<target .and. yelm(6)<target) &
    iMPIcut_eta(1,ispec)=.true.

! right cut-plane in the current slice along Y = constant (Ymax of this slice)
! and add geometrical tolerance

  target = (iproc_eta+1)*sizeslice - TOLERANCE_DEGREES
  if(yelm(3)>target .and. yelm(4)>target .and. yelm(7)>target .and. yelm(8)>target) &
    iMPIcut_eta(2,ispec)=.true.

  end subroutine get_flags_boundaries_cubed

