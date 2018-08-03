!======================================================================================
!
!  CUBED_SPHERE projection from the global SEM3D of Dimitri Komatitsch and Jeroen Tromp
!
!======================================================================================
subroutine cubed_geo(xelm,yelm,zelm,ANGULAR_WIDTH_XI_IN_DEGREES, &
           ANGULAR_WIDTH_ETA_IN_DEGREES,rotation_matrix)

  use constants
  implicit none
!  include "constants.h"

  double precision xelm(NGNOD_EIGHT_CORNERS)
  double precision yelm(NGNOD_EIGHT_CORNERS)
  double precision zelm(NGNOD_EIGHT_CORNERS)

  double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES
  ! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix
  double precision, parameter :: degrad=PI/180.d0
  !WENBO
  !double precision, parameter :: r_earth=6371000.0
  double precision ::x,y,x_sphere,y_sphere,z_sphere
  double precision ::r,rt
  double precision, dimension(NDIM) :: vector_ori,vector_rotated
  double precision ::gamma

  integer ::ia,i,j
  ! loop over the NGNOD_EIGHT_CORNERS nodes
  do ia=1,NGNOD_EIGHT_CORNERS 
     r=R_TOP_BOUND+zelm(ia)
     x = dtan((xelm(ia)-ANGULAR_WIDTH_XI_IN_DEGREES/2.d0)*degrad)
     y = dtan((yelm(ia)-ANGULAR_WIDTH_ETA_IN_DEGREES/2.d0)*degrad)
     gamma = ONE / dsqrt(ONE + x*x + y*y)
     rt=r*gamma
     x_sphere = -y*rt
     y_sphere = x*rt
     z_sphere = rt

!WENBO
!to be deleted
!     if(dsqrt(x_sphere**2+y_sphere**2+z_sphere**2)<6320999.0)&
!      print *,'Error',r,R_TOP_BOUND,zelm(ia),dsqrt(x_sphere**2+y_sphere**2+z_sphere**2)
 
   ! rotate
     vector_ori(1) = x_sphere
     vector_ori(2) = y_sphere
     vector_ori(3) = z_sphere
     do i = 1,NDIM
        vector_rotated(i) = ZERO
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     xelm(ia) = vector_rotated(1)
     yelm(ia) = vector_rotated(2)
     zelm(ia) = vector_rotated(3)
   end do
end subroutine cubed_geo
