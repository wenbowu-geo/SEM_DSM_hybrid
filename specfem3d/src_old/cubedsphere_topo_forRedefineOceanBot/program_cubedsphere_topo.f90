!=================================================================
!                       Cubedsphere_topography
!           --------------------------------
!
!                   author:Wenbo Wu
!                          USTC
!                       April 2013
!
!   This program is to convert the topography in the lat-lon coordinate 
!   system into that in the cubed sphere coordinate.

!=================================================================
!########The code converts the real topography in longitude-latitude system into
!that
!in a cubed spehre xi-eta system. Then the file
!"your_example_folder/DATA/meshfem3D_files/real_bathymetry_topography" is used
!in the step xmeshfem3D, in which the elements are classified as elastic or
!acoustic types.


 ! reads in arguments

 !The first input argument arge(1) is the name of topography model, as a
 !function
 !of longitude and latitude coordinates. This file has the below format

 ! nlongitude nlatitude
 ! origin_longitude origin_latitude
 ! space_longitude space_latitude
!  topography(lat0,lon0) topography(lat0,lon0+dlon) topography(lat0,lon0+2*dlon)
!  topography(lat0,lon0+3*dlon) ...
!  topography(lat0+dlat,lon0) topography(lat0+dlat,lon0+dlon) ...
!  ...
!  topography(lat0+(nlongitude-1)*dlat,lon0) ...
!  topography(lat0+(nlongitude-1)*dlat,lon0+(nlon-1)*dlon) 





 !The remained  input arguments are
 !     arge(2) - nxi
 !     arge(3) - neta
 !     arge(4) - origin_xi
 !     arge(5) - origin_eta
 !     arge(6) - space_xi
 !     arge(6) - space_eta


 !USAGE: ../../bin/xcubedsphere_topo_forRedineOcean
 !./DATA/meshfem3D_files/SA_real_topo.dat 540 540 -0.2        -0.2d0      0.01
 !0.01



 !The above command converts the topography "SA_real_topo.dat" model in the
 !lon-lat
 !coordinate system into that in a cubed sphere (xi-eta) coordinate system.
 !The converted topography are saved in the file 
 !your_example_folder/DATA/meshfem3D_files/real_bathymetry_topography. This file
 !name is used
 !to determine whether the element is solid or fluid (ocean).


!WENBOWU
!Data: 2018-02-05

  program  xcubedsphere_topo

!  call initialize_cubed_topo()
  call cubedsphere_topo()

  end program xcubedsphere_topo
