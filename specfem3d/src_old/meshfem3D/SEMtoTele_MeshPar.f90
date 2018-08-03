!=====================================================================
!

!module constants

!  include "constants.h"

!end module constants

!=====================================================================
module SEMtoTele_MeshPar
use constants

   double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
                   CENTER_LATITUDE_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH
   integer::NX_TOPOTAPER,NY_TOPOTAPER,NX_NOTOPO,NY_NOTOPO
   integer::TELE_IXLOW,TELE_IXHIGH,TELE_IYLOW,TELE_IYHIGH,TELE_IRTOP,TELE_IRBOT
   integer, parameter ::Imain_SEMtoTele=1011


   integer ::TeleEle_nxLow,TeleEle_nxHigh,TeleEle_nyLow,&
             TeleEle_nyHigh,TeleEle_nrTop,TeleEle_nrBot
   integer,dimension(:),allocatable  ::TeleEle_xLow,TeleEle_xLowReg
   integer,dimension(:),allocatable  ::TeleEle_xHigh,TeleEle_xHighReg
   integer,dimension(:),allocatable  ::TeleEle_yLow,TeleEle_yLowReg
   integer,dimension(:),allocatable  ::TeleEle_yHigh,TeleEle_yHighReg
   integer,dimension(:),allocatable  ::TeleEle_rTop,TeleEle_rTopReg
   integer,dimension(:),allocatable  ::TeleEle_rBot,TeleEle_rBotReg

   double precision, dimension(NDIM,NDIM) :: rotation_matrix

   double precision, dimension(:,:,:,:), allocatable :: xstore_cubedsph,ystore_cubedsph,zstore_cubedsph
   double precision, dimension(:,:), allocatable :: nodes_coords_cubedsph,nodes_coords_cubedsph_old


!   logical ::LOW_RESOLUTION
!   integer ::NSTEP_BETWEEN_OUTPUTBOUND
!   integer ::DECIMATE_COUPLING
!   integer ::NPOINTS_PER_PACK
end module SEMtoTele_MeshPar


