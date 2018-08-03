
!########The code converts the real topography in longitude-latitude system into that
!in a cubed spehre xi-eta system. Then the file
!"your_example_folder/DATA/meshfem3D_files/real_bathymetry_topography" is used
!in the step xmeshfem3D, in which the elements are classified as elastic or acoustic types.


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
  double precision:: orig_x_interface,orig_y_interface
  double precision :: orig_lat_interface,orig_lon_interface
  double precision :: spacing_x_interface,spacing_y_interface
  double precision:: spacing_lat_interface,spacing_lon_interface
  character(len=70) INTERFACES_FILE,CAVITY_FILE
  character(len=70), parameter::interface_file="real_bathymetry_topography"
  character(len=70)::topo_latlon_interface_file
  integer i_interface ! ipoint_current
  integer number_of_interfaces
  integer ::latlon_number_of_interfaces
  integer :: nx_interface,ny_interface
  integer nlat_interface,nlon_interface
  double precision:: topography_basement

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

 topo_latlon_interface_file=arg(1)
 read(arg(2),*) nx_interface
 read(arg(3),*) ny_interface
 read(arg(4),*) orig_x_interface
 read(arg(5),*) orig_y_interface
 read(arg(6),*) spacing_x_interface
 read(arg(7),*) spacing_y_interface



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
!            call read_value_double_precision(45,DONT_IGNORE_JUNK,topo_latlon_interface(ilat,ilon),'Z_INTERFACE_TOP')
!           read(45,*) topo_latlon_interface(ilat,ilon)
        enddo
    enddo

! goes from LatMin to LatMax
    do ilat=1,nlat_interface
       read(45,*) topo_latlon_interface(ilat,:)
    enddo

! goes from LatMax to LatMin
    do ilat=nlat_interface,1,-1
!       read(45,*) topo_latlon_interface(ilat,:)
    enddo

    
    close(45)
!    call smooth_interface(topo_latlon_interface,nlat_interface(i_interface),nlon_interface(i_interface))



!******************************************************************************
!**************Convert topography(lon,lat) into topography(xi,eta) ************
    allocate(xy_interface(nx_interface,ny_interface))
    !print *,'nx ny',nx_interface,ny_interface
    do ix=1,nx_interface
      do iy=1,ny_interface
        x_current=dble(ix-1)*spacing_x_interface+orig_x_interface
        y_current=dble(iy-1)*spacing_y_interface+orig_y_interface

        call cubedspheretoxyz(x_current,y_current,x,y,z,rotation_matrix,ANGULAR_WIDTH_XI_IN_DEGREES, &
                              ANGULAR_WIDTH_ETA_IN_DEGREES)
        call xyztolatlon(x,y,z,lat,lon)
        call search_topo(lat,lon,xy_interface(ix,iy),topo_latlon_interface,nlat_interface,nlon_interface,&
                         orig_lat_interface,orig_lon_interface,spacing_lat_interface,spacing_lon_interface)



!***********************No taper here, because this file is used to determine**********
!***********************acoustic(ocean) or elastic(land)**********************



!************************add taper******************************************
!        if(x_current.le.belowx_notopo.or.y_current.le.belowy_notopo.or.&
!           x_current.ge.abovex_notopo.or.y_current.ge.abovey_notopo) then
!            xy_interface(ix,iy)=topography_basement
!        else
!           if(x_current.le.belowx_taper) &
!             xy_interface(ix,iy)=topography_basement+ &
!                (xy_interface(ix,iy)-topography_basement)*(x_current-belowx_notopo)/xtaper_width
!           if(x_current.ge.abovex_taper) &
!             xy_interface(ix,iy)=topography_basement+ &
!                (xy_interface(ix,iy)-topography_basement)*(abovex_notopo-x_current)/xtaper_width
!           if(y_current.le.belowy_taper) &
!             xy_interface(ix,iy)=topography_basement+ &
!                (xy_interface(ix,iy)-topography_basement)*(y_current-belowy_notopo)/ytaper_width
!           if(y_current.ge.abovey_taper) &
!             xy_interface(ix,iy)=topography_basement+ &
!                (xy_interface(ix,iy)-topography_basement)*(abovey_notopo-y_current)/ytaper_width
!        end if

      end do !loop for iy
    end do  !loop for ix

    call save_topo(nx_interface,ny_interface,xy_interface, &
            interface_file,orig_x_interface,orig_y_interface,&
            spacing_x_interface,spacing_y_interface)




    deallocate(xy_interface)
    deallocate(topo_latlon_interface)


  end subroutine cubedsphere_topo

