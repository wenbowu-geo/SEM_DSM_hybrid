subroutine cubedspheretoxyz(x_cubed,y_cubed,x,y,z,rotation_matrix,ANGULAR_WIDTH_XI_IN_DEGREES,&
                            ANGULAR_WIDTH_ETA_IN_DEGREES)
  implicit none
  include 'constants.h'


  double precision, intent(in):: ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES
  double precision, intent(in), dimension(NDIM,NDIM) :: rotation_matrix
  double precision, intent(in) ::x_cubed,y_cubed
  double precision, intent(out) ::x,y,z

  !auxiliary variables
  double precision, parameter :: degrad=PI/180.d0
  double precision, parameter :: raddeg=180.d0/PI
  !double precision, parameter :: r_earth=6371000.00
  integer ::i,j
  double precision ::x_cubedsphere,y_cubedsphere
  double precision ::r,rt,gamma
  double precision, dimension(NDIM) :: vector_ori,vector_rotated


  x_cubedsphere=tan((x_cubed-ANGULAR_WIDTH_XI_IN_DEGREES/2.d0)*degrad)
  y_cubedsphere=tan((y_cubed-ANGULAR_WIDTH_ETA_IN_DEGREES/2.d0)*degrad)
  gamma = ONE / sqrt(ONE + x_cubedsphere*x_cubedsphere + y_cubedsphere*y_cubedsphere)
  rt=R_TOP_BOUND*gamma
  vector_ori(1)=-y_cubedsphere*rt
  vector_ori(2)=x_cubedsphere*rt
  vector_ori(3)=rt

  !rotate
  do i = 1,NDIM 
      vector_rotated(i) = ZERO
      do j = 1,NDIM
          vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
      enddo
  enddo
  x=vector_rotated(1)
  y=vector_rotated(2)
  z=vector_rotated(3)
end subroutine cubedspheretoxyz




subroutine xyztolatlon(x,y,z,lat,lon)
  implicit none
  include 'constants.h'

  double precision, intent(in) ::x,y,z
  double precision, intent(out) ::lat,lon
  
  

!auxiliary variables
  double precision ::theta,phi
  double precision, parameter :: raddeg=180.d0/PI
!  double precision, parameter :: r_earth=6371000.00


  theta=acos(z/R_TOP_BOUND)
  phi=acos(x/(R_TOP_BOUND*dsin(theta))*(1-1.e-7))
  theta=theta*raddeg
  phi=phi*raddeg
!0-360deg
  if(y<-1.e-14) phi=-phi+360.0
!-180 to 180deg
!  if(y<-1.e-14) phi=-phi


  lat=-theta+90.0
  lon=phi
end subroutine xyztolatlon
