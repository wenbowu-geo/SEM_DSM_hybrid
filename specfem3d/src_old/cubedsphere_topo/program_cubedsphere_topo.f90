!=================================================================
!                       Cubedsphere_topo
!           --------------------------------
!
!                   author:Wenbo Wu
!                          USTC
!                       April 2013
!
!   This program is to convert the topography in the lat-lon coordinate 
!   system into that in the cubed-sphere sphere coordinate.

!=================================================================
!########The code use linear interpolation to convert the topography in the
!longitude-latitude system into that
!in a cubed spehre xi-eta system
!WENBOWU
!Data: 2018-02-05


!********************************Input
!argument**********************************
 !The first input argument arg(1) specifies which individual interface we are
 !working on. All the optional interfeaces are listed in the file
 !your_example_folder/DATA/meshfem3D_files/interfaces.dat. Usually, the last one
 !is free
 !surface and the last but one could be sediment base or ocean bottom.

 !The second input argument arge(2) is the name of topography model a function
 !of longitude and latitude coordinates.

 !The third  input argument arge(3) is the basement value of this interface.

 !The fouth input argument arge(4) is the minimum allowed topography for this 
 !particular interface. 
 !The fouth input argument arge(5) is the maximum allowed topography for this 
 !particular interface.
 !For some cases, specifying an minimum and/or maximum allowed topography is
 !useful. 
 !For example, you may need to create an ocean layer in the simualtion and the 
 ! topography/bathymetry is separated to two interfaces. The top one 
 !is (land) topography and 0 meter (a minimum allowed topography) in the ocean
 !area. 
 !The other interface is ocean bottom and artifical interface (maximum allowed
 !topography) 
 !beneath land. By doing this and choosing a proper truncated topography or
 !bathymetry, 
 !we create a space for the ocean layer. Otherwise, this layer would be too thin
 !or the two interfaces would intersect with each other.
!********************************************************************************




!****************************************Example************************************
 ! If the file your_example_folder/DATA/meshfem3D_files/interfaces.dat is as
 ! below

   !  # number of interfaces
   !  3
   !#
   !# We describe each interface below, structured as a 2D-grid, with several
   !parameters :
   !# number of points along XI and ETA, minimal XI ETA coordinates
   !# and spacing between points which must be constant.
   !# Then the records contain the Z coordinates of the NXI x NETA points.
   !#
   !# interface number 1 (topography, top of the mesh)
   !# SUPPRESS_UTM_PROJECTION  NXI  NETA LONG_MIN   LAT_MIN    SPACING_XI
   !SPACING_ETA

   !.true.                   110 110 -0.2        -0.2      0.05d0    0.05d0
   ! topo_Moho.dat
   !.true.                   110 110 -0.2        -0.2      0.05d0    0.05d0
   !topo_Midcrust.dat
   !.true.                   320 320 -0.1        -0.1      0.01d0    0.01d0
   !topo_free_surf.dat
   !7
   !2
   !2


!~~~~~~~Free surface topography~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! go to your work folder and run the command
 ! ../../bin/xtopo_lonlat_to_cubedsphere 3
 ! ./DATA/meshfem3D_files/ETOPO1_topography 0.0 -12000.0 10000.0

 !The above command converts the topography "ETOPO1_topography" in lon-lat
 !coordinate system into that in a cubed sphere (xi-eta) coordinate system.
 !The converted topography are saved in the file
 !your_example_folder/DATA/meshfem3D_files/topo_free_surf.dat. This file name is
 !the third
 !interface given in the file
 !your_example_folder/DATA/meshfem3D_files/interfaces.dat
 !(See above). The basement of free surface is 0.0, so the topography toward the
 !boundary of SEM box is tapered to 0.0.

! This file ETOPO1_topography has the below format
!-------------- the format of ETOPO1_topography----------------------
!  nlongitude nlatitude
!  origin_longitude origin_latitude
!  space_longitude space_latitude
!  topography(lat0,lon0) topography(lat0,lon0+dlon) topography(lat0,lon0+2*dlon)
!  topography(lat0,lon0+3*dlon) ...
!  topography(lat0+dlat,lon0) topography(lat0+dlat,lon0+dlon) ...
!  ...
!  topography(lat0+(nlongitude-1)*dlat,lon0) ... topography(lat0+(nlongitude-1)*dlat,lon0+(nlon-1)*dlon)

!------------------------------------------------------------------
! Be careful!!! (1) The first three lines privide some basic information of this
! file.
! After that, each line is compoded of topography at one latitude and all the
! longitudes.
! Or you can change the lines 456-463 script to fit the format of your model.
! (2)
! the longitude should be between 0-360deg.

! The minimum and maximum allowed topogrpy/bathymetry are -12000.0 (-12km) and
! 10000.0 (10 km) respectively.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~Moho discontinuity~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! go to your work folder and run the command
 !  ../../bin/xtopo_lonlat_to_cubedsphere 1
 !  ./DATA/meshfem3D_files/Moho_topography_latlon -35000.0 -60000.0 -10000.0

 !The above command converts the topography "Moho_topography_latlon" in lon-lat
 !coordinate system into that in a cubed sphere (xi-eta) coordinate system.
 !The converted topography are saved in the file
 !your_example_folder/DATA/meshfem3D_files/topo_Moho.dat. This file name is the
 !first
 !interface given in the file
 !your_example_folder/DATA/meshfem3D_files/interfaces.dat
 !The basement of Moho is -35000.0 (35 km, it is the crust thickness in iasp91
 !model),
 !so the Moho topography toward the boundary of SEM box is tapered to -35000.0.
 !The minimum and maximum allowed Moho topogrphy is -60000.0 (-60 km) and
 !-10000.0 (-10 km).


 !Keep in mind, the basement of discontinuity must be consistent with the
 !background velocity model "your_example_folder/DATA/dsm_model_input". This is
 !reauired by the coupling theory. However, the ocean is an exception, because
 !you may want to exclude
 !ocean in the dsm_imput_file, but include it in the SEM model.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  program  xcubedsphere_topo

!  call initialize_cubedsphere_topo()
  call cubedsphere_topo()

  end program xcubedsphere_topo
