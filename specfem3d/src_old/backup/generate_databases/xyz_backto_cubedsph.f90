!======================================================================================
!
!  Refer to - CUBED_SPHERE projection from the global SEM3D of Dimitri Komatitsch and Jeroen Tromp
!
!======================================================================================
subroutine xyz_backto_cubedsph(x,y,z,x_cubedsph,y_cubedsph,z_cubedsph,&
           ANGULAR_WIDTH_XI_IN_DEGREES, &
           ANGULAR_WIDTH_ETA_IN_DEGREES,rotation_matrix_back_cubedsph)

  use constants
  implicit none
!  include "constants.h"

  double precision,intent(in):: x,y,z
  double precision,intent(out):: x_cubedsph,y_cubedsph,z_cubedsph

  double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES
  ! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix_back_cubedsph
  double precision, parameter :: degrad=PI/180.d0
  !WENBO
  !double precision, parameter :: r_earth=6371000.0
  double precision, dimension(NDIM) :: vector_ori,vector_rotated

  integer ::ia,i,j
  vector_ori(1) = x
  vector_ori(2) = y
  vector_ori(3) = z
  do i = 1,NDIM
     vector_rotated(i) = ZERO
     do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix_back_cubedsph(i,j)*vector_ori(j)
     enddo
  enddo
!z coordinates is the raius rather than the real z in cubed sphere coordinates
!system!
!  z_cubedsph=vector_rotated(3)
  z_cubedsph=dsqrt(x**2+y**2+z**2)
!  if(z_cubedsph<6320999.0) print *,"Error",z_cubedsph,x,y,z
  y_cubedsph=atan(-vector_rotated(1)/vector_rotated(3))
  x_cubedsph=atan(vector_rotated(2)/vector_rotated(3))

!make y_cubedsph and x_cubedsph in unit of deg, rather than rad. Because in the
!subrountine pml_set_local_dampingcoeff.f90, 

end subroutine xyz_backto_cubedsph
