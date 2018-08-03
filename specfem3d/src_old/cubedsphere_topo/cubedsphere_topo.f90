
!########The code use linear interpolation to convert the topography in the longitude-latitude system into that
!in a cubed spehre xi-eta system
!WENBOWU
!Data: 2018-02-05


!********************************Input argument**********************************
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
 !For some cases, specifying an minimum and/or maximum allowed topography is useful. 
 !For example, you may need to create an ocean layer in the simualtion and the 
 ! topography/bathymetry is separated to two interfaces. The top one 
 !is (land) topography and 0 meter (a minimum allowed topography) in the ocean area. 
 !The other interface is ocean bottom and artifical interface (maximum allowed topography) 
 !beneath land. By doing this and choosing a proper truncated topography or bathymetry, 
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
 ! ../../bin/xtopo_lonlat_to_cubedsphere 3 ./DATA/meshfem3D_files/ETOPO1_topography 0.0 -12000.0 10000.0

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
 !The minimum and maximum allowed Moho topogrphy is -60000.0 (-60 km) and -10000.0 (-10 km).


 !Keep in mind, the basement of discontinuity must be consistent with the
 !background velocity model "your_example_folder/DATA/dsm_model_input". This is
 !reauired by the coupling theory. However, the ocean is an exception, because
 !you may want to exclude
 !ocean in the dsm_imput_file, but include it in the SEM model.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





  subroutine cubedsphere_topo
  implicit none
  include "constants.h"
!==========================================================



! auxiliary variables to generate the mesh
  integer ix,iy
  double precision x_current,y_current !,z_top,z_bot
  double precision x,y,z
  integer ilat,ilon
  double precision lat,lon







! parameters read from parameter file
  integer NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE


  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX
  double precision Z_DEPTH_BLOCK !,Z_BASEMENT_SURFACE,Z_DEPTH_MOHO
  double precision LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX
  double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
                   CENTER_LATITUDE_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH
  ! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix

  !logical TOPOGRAPHY,MOHO_MAP_LUPEI
  !logical BASEMENT_MAP,HAUKSSON_REGIONAL_MODEL,HARVARD_3D_GOCAD_MODEL,IMPOSE_MINIMUM_VP_GOCAD

  logical SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH
  logical CUBED_SPHERE_PROJECTION

! Mesh files for visualization
  logical CREATE_ABAQUS_FILES,CREATE_DX_FILES

! doublings parameters
  integer NDOUBLINGS
  integer, dimension(2) :: ner_doublings

  character(len=256) OUTPUT_FILES,LOCAL_PATH !,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA

! this for all the regions
  integer NSPEC_AB,NGLOB_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX

  double precision min_elevation,max_elevation
  !double precision min_elevation_all,max_elevation_all


  !parameter for SEMtoTele
    integer ::nx_TopoTaper,ny_TopoTaper,nx_notopo,ny_notopo
    integer ::Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,Tele_irTop,Tele_irBot
    integer ::TeleEle_ixLow,TeleEle_ixHigh,TeleEle_iyLow,TeleEle_iyHigh,TeleEle_irTop,TeleEle_irBot
    integer ::TeleEle_nxLow,TeleEle_nxHigh,TeleEle_nyLow,TeleEle_nyHigh,TeleEle_nrTop,TeleEle_nrBot
   
    double precision::belowx_notopo,belowx_taper,abovex_notopo,abovex_taper
    double precision::belowy_notopo,belowy_taper,abovey_notopo,abovey_taper
    double precision::xtaper_width,ytaper_width


! for tapered basement map
  !integer icorner_x,icorner_y
  !integer iz_basement
  !double precision x_corner,y_corner
  !character(len=256) BASEMENT_MAP_FILE

! interfaces parameters
  double precision,dimension(:,:),pointer:: xy_interface
  double precision,dimension(:,:),pointer:: topo_latlon_interface
  double precision,dimension(:),pointer:: orig_x_interface,orig_y_interface
  double precision :: orig_lat_interface,orig_lon_interface
  double precision,dimension(:),pointer:: spacing_x_interface,spacing_y_interface
  double precision:: spacing_lat_interface,spacing_lon_interface
  character(len=70) INTERFACES_FILE,CAVITY_FILE
  character(len=70), dimension(:), pointer::interface_file
  character(len=70)::topo_latlon_interface_file
  integer i_interface ! ipoint_current
  integer number_of_interfaces
  integer ::latlon_number_of_interfaces
  integer, dimension(:),pointer :: nx_interface,ny_interface
  integer nlat_interface,nlon_interface
  integer ::ith_interface_workon
  double precision:: topography_basement
  double precision:: max_topography_truncation,min_topography_truncation

! subregions parameters
  integer NSUBREGIONS
!  definition of the different regions of the model in the mesh (nx,ny,nz)
!  #1 #2 : nx_begining,nx_end
!  #3 #4 : ny_begining,ny_end
!  #5 #6 : nz_begining,nz_end
!     #7 : material number
  integer, dimension(:,:), pointer :: subregions

! material properties
  integer NMATERIALS
! first dimension  : material_id
! second dimension : #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id
  double precision , dimension(:,:), pointer :: material_properties


!arguments 
  integer ::iarg
  character(len=MAX_STRING_LEN) :: arg(9)



! ************** PROGRAM STARTS HERE **************
 do iarg = 1, command_argument_count()
    call get_command_argument(iarg,arg(iarg))
 enddo

 read(arg(1),*) ith_interface_workon
 topo_latlon_interface_file=arg(2)
 read(arg(3),*) topography_basement
 read(arg(4),*) min_topography_truncation
 read(arg(5),*) max_topography_truncation


 if(topography_basement>max_topography_truncation) then
   print *,"Error, check your inputs, min_topography_truncation>max_topography_truncation!!"
   call exit()
 end if

 if(topography_basement>max_topography_truncation) then
   print *, 'The topography is truncated below',max_topography_truncation,&
        ' meters, but the topography_basement is even higher than the that.'
   print *,"Error, check your input, topography_basement>max_topography_truncation!!"
   call exit()
 end if

 if(topography_basement<min_topography_truncation) then
   print *, 'The topography is truncated above', min_topography_truncation,&
         ' meters, but the topography_basement is even lower than the that.'
   print *,"Error, check your input, topography_basement>max_topography_truncation!!"
   call exit()
 end if


! get the base pathname for output files
!  call get_value_string(OUTPUT_FILES_BASE, 'OUTPUT_FILES', OUTPUT_FILES_BASE(1:len_trim(OUTPUT_FILES_BASE)))

 open(unit=IMAIN,file=OUTPUT_FILES_BASE(1:len_trim(OUTPUT_FILES_BASE))//'/output_topo_convert.txt',status='unknown')

 write(IMAIN,*)
 write(IMAIN,*) '******************************************'
 write(IMAIN,*) '*** Specfem3D MPI Topography_convert - f90 version ***'
 write(IMAIN,*) '******************************************'
 write(IMAIN,*)



  call read_parameter_file(LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK,CUBED_SPHERE_PROJECTION,&
        NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE,LOCAL_PATH,SUPPRESS_UTM_PROJECTION,&
        INTERFACES_FILE,CAVITY_FILE,NSUBREGIONS,&
        USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)



  if(.not.CUBED_SPHERE_PROJECTION) stop 'This program is only for the CUBED SPHERE PROJECTION=.TRUE.'


  write(IMAIN,*)
  write(IMAIN,*) 'Reading the taper information from file ', &
       MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES)) &
       //'SEMtoTele_Par_file'

  call read_parameter_SEMtoTele(IMAIN,nx_TopoTaper,ny_TopoTaper,&
         nx_notopo,ny_notopo,Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,&
         Tele_irTop,Tele_irBot,NEX_XI,NEX_ETA,10000000,NDOUBLINGS,USE_REGULAR_MESH,&
         ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
         CENTER_LATITUDE_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,& 
         rotation_matrix,ner_doublings,0)
!**********************notopography and taper range*************************
         xtaper_width=nx_TopoTaper*ANGULAR_WIDTH_XI_IN_DEGREES/NEX_XI
         belowx_notopo=nx_notopo*ANGULAR_WIDTH_XI_IN_DEGREES/NEX_XI
         belowx_taper= belowx_notopo+xtaper_width
         abovex_notopo=ANGULAR_WIDTH_XI_IN_DEGREES*(1.0-dble(nx_notopo)/NEX_XI)
         abovex_taper=abovex_notopo-xtaper_width

         ytaper_width=ny_TopoTaper*ANGULAR_WIDTH_ETA_IN_DEGREES/NEX_ETA
         belowy_notopo=ny_notopo*ANGULAR_WIDTH_ETA_IN_DEGREES/NEX_ETA
         belowy_taper=belowy_notopo+ytaper_width
         abovey_notopo=ANGULAR_WIDTH_ETA_IN_DEGREES*(1.0-dble(ny_notopo)/NEX_ETA)
         abovey_taper=abovey_notopo-ytaper_width
  write(IMAIN,*)
  write(IMAIN,*) 'Taper information'
  write(IMAIN,*) 'Below x=',belowx_notopo,'and above ',abovex_notopo,' topography is zero.'
  write(IMAIN,*) 'For x:',belowx_notopo,'-',belowx_taper,'deg and',abovex_taper,'-',abovex_notopo,&
                 'deg is topography tapering range.'
  write(IMAIN,*) 'Below y=',belowy_notopo,'and above ',abovey_notopo,' topography is zero.'
  write(IMAIN,*) 'For y',belowy_notopo,'-',belowy_taper,'deg and',abovey_taper,'-',abovey_notopo,&
                 'deg is topography tapering range.'


!*****************************************************************************


! get interface data from external file to count the spectral elements along Z
  write(IMAIN,*) 'Reading xy_interface data from file ',&
        MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES))//INTERFACES_FILE(1:len_trim(INTERFACES_FILE)), &
        ' to count the spectral elements'

  open(unit=IIN,file=MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES))//INTERFACES_FILE,status='old')


! read number of interfaces
  call read_value_integer(IIN,DONT_IGNORE_JUNK,number_of_interfaces,'NINTERFACES')
  if(number_of_interfaces < 1) stop 'not enough interfaces (minimum is 1, for topography)'

  allocate(interface_file(number_of_interfaces))
  allocate(nx_interface(number_of_interfaces))
  allocate(ny_interface(number_of_interfaces))
  allocate(orig_x_interface(number_of_interfaces))
  allocate(orig_y_interface(number_of_interfaces))
  allocate(spacing_x_interface(number_of_interfaces))
  allocate(spacing_y_interface(number_of_interfaces))


! loop on all the interfaces
  do i_interface = 1,number_of_interfaces

     call read_interface_parameters(IIN,SUPPRESS_UTM_PROJECTION,interface_file(i_interface), &
          nx_interface(i_interface),ny_interface(i_interface),orig_x_interface(i_interface),orig_y_interface(i_interface),&
          spacing_x_interface(i_interface),spacing_y_interface(i_interface))
     write(IMAIN,*) 'the ',i_interface,'interface: nx_interface=',nx_interface(i_interface),'ny_interface=',ny_interface(i_interface)

     if((nx_interface(i_interface) < 2) .or.(ny_interface(i_interface) < 2)) stop 'not enough interface points (minimum is 2x2)'
  enddo


  close(IIN)



  min_elevation = +HUGEVAL
  max_elevation = -HUGEVAL


!*********************************************************************************
!*******************read the input topography(longitude,latitude) model************
!  write(IMAIN,*)
!  write(IMAIN,*) 'Reading topo_latlon_interface data from file ', &
!        MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES)) &
!        //INTERFACES_FILE(1:len_trim(INTERFACES_FILE))
!  write(IMAIN,*)

!  open(unit=IIN,file=MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES)) &
!          //"latlon_"//INTERFACES_FILE,status='old')


  ! read number of interfaces
!  call read_value_integer(IIN,DONT_IGNORE_JUNK,latlon_number_of_interfaces,'NINTERFACES')


!    call read_interface_parameters(IIN,SUPPRESS_UTM_PROJECTION,topo_latlon_interface_file,&
!             nlon_interface,nlat_interface,&
!             orig_lon_interface,orig_lat_interface,spacing_lon_interface,spacing_lat_interface)
!    allocate(topo_latlon_interface(nlat_interface,nlon_interface))


    ! loop on all the points describing this interface
    open(unit=45,file=topo_latlon_interface_file,status='old')
    read(45,*) nlon_interface,nlat_interface

    allocate(topo_latlon_interface(nlat_interface,nlon_interface))
    read(45,*) orig_lon_interface,orig_lat_interface
    read(45,*) spacing_lon_interface,spacing_lat_interface

!The most inner loop for lat and next for lon
    do ilon=1,nlon_interface
        do ilat=1,nlat_interface
!            call read_value_double_precision(45,DONT_IGNORE_JUNK,topo_latlon_interface(ilat,ilon),'Z_INTERFACE_TOP')
!           read(45,*) topo_latlon_interface(ilat,ilon)
        enddo
    enddo

!The most inner loop for lat and next for lon
    do ilat=1,nlat_interface
        do ilon=1,nlon_interface
!            call
!            read_value_double_precision(45,DONT_IGNORE_JUNK,topo_latlon_interface(ilat,ilon),'Z_INTERFACE_TOP')
!           read(45,*) topo_latlon_interface(ilat,ilon)
        enddo
    enddo

! goes from LatMin to LatMax
    do ilat=1,nlat_interface
       read(45,*) topo_latlon_interface(ilat,:)
    enddo

! goes from LatMax to LatMin
!    do ilat=nlat_interface,1,-1
!       read(45,*) topo_latlon_interface(ilat,:)
!    enddo


    close(45)
!    call smooth_interface(topo_latlon_interface,nlat_interface(i_interface),nlon_interface(i_interface))



!******************************************************************************
!**************Convert topography(lon,lat) into topography(xi,eta) ************
  do i_interface=ith_interface_workon,ith_interface_workon  !1,number_of_interfaces
    allocate(xy_interface(nx_interface(i_interface),ny_interface(i_interface)))
    !print *,'nx ny',nx_interface(i_interface),ny_interface(i_interface)
    do ix=1,nx_interface(i_interface)
      do iy=1,ny_interface(i_interface)
        x_current=dble(ix-1)*spacing_x_interface(i_interface)+orig_x_interface(i_interface)
        y_current=dble(iy-1)*spacing_y_interface(i_interface)+orig_y_interface(i_interface)

        call cubedspheretoxyz(x_current,y_current,x,y,z,rotation_matrix,ANGULAR_WIDTH_XI_IN_DEGREES, &
                              ANGULAR_WIDTH_ETA_IN_DEGREES)
        call xyztolatlon(x,y,z,lat,lon)
        call search_topo(lat,lon,xy_interface(ix,iy),topo_latlon_interface,nlat_interface,nlon_interface,&
                         orig_lat_interface,orig_lon_interface,spacing_lat_interface,spacing_lon_interface)
!************************add taper*******************************
        if(x_current.le.belowx_notopo.or.y_current.le.belowy_notopo.or.&
           x_current.ge.abovex_notopo.or.y_current.ge.abovey_notopo) then
            xy_interface(ix,iy)=topography_basement
        else
           if(x_current.le.belowx_taper) &
             xy_interface(ix,iy)=topography_basement+ &
                (xy_interface(ix,iy)-topography_basement)*(x_current-belowx_notopo)/xtaper_width
           if(x_current.ge.abovex_taper) &
             xy_interface(ix,iy)=topography_basement+ &
                (xy_interface(ix,iy)-topography_basement)*(abovex_notopo-x_current)/xtaper_width
           if(y_current.le.belowy_taper) &
             xy_interface(ix,iy)=topography_basement+ &
                (xy_interface(ix,iy)-topography_basement)*(y_current-belowy_notopo)/ytaper_width
           if(y_current.ge.abovey_taper) &
             xy_interface(ix,iy)=topography_basement+ &
                (xy_interface(ix,iy)-topography_basement)*(abovey_notopo-y_current)/ytaper_width
        end if

!*******************************Truncating topography higher than max_topography_truncation******
        if(xy_interface(ix,iy)>max_topography_truncation) xy_interface(ix,iy)=max_topography_truncation
        if(xy_interface(ix,iy)<min_topography_truncation) xy_interface(ix,iy)=min_topography_truncation


      end do !loop for iy
    end do  !loop for ix

    call save_topo(nx_interface(i_interface),ny_interface(i_interface),xy_interface,interface_file(i_interface))
    deallocate(xy_interface)
    deallocate(topo_latlon_interface)
  end do !loop for i_interface


  end subroutine cubedsphere_topo

