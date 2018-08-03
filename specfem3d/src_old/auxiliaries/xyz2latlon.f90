subroutine  xyz2latlon(x,y,z,lat,lon)
implicit none
include "constants.h"
double precision::x,y,z
double precision::lat,lon,Radi,theta,phi
!double precision::PI= 3.1415926535897932d0

Radi=dsqrt(x**2+y**2+z**2)
  theta=dacos(z/Radi)
  if(y*dsin(Theta)<0.0) then
       Phi=2*PI-dacos(x*(1.d0-1.e-10)/(Radi*dsin(Theta)))
  else
       Phi=dacos(x*(1.d0-1.e-10)/(Radi*dsin(Theta)))
  end if 
  lon=Phi/PI*180.0
!  if(ELLIPTICITY) then
!     dcost = dcos(ReceiverInfo(i)%theta)
!     p20 = 0.5d0*(3.0d0*dcost*dcost-1.0d0)
!     call spline_evaluation(rspl,espl,espl2,nspl,R_EARTH,ell)
!     Rearth  =  R_EARTH*(1.0d0-(2.0d0/3.0d0)*ell*p20)
!     lat= 180.0*PI*atan(tan(PI/2-theta)/0.99329534d0)
!  else  
!     Rearth  =  R_EARTH
!     lat= 180/PI*(PI/2-theta)
!  end if

  lat= 180.0/PI*(PI/2-theta)

end subroutine
