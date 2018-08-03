subroutine search_topo(lat,lon,topography,latlon_interface,nlat_interface,nlon_interface,&
orig_lat_interface,orig_lon_interface,spacing_lat_interface,spacing_lon_interface)
  implicit none
  include 'constants.h'
  
  double precision,intent(in):: lat,lon
  double precision, intent(out):: topography
  integer, intent(in)::nlat_interface,nlon_interface
  double precision, intent(in)::orig_lat_interface,orig_lon_interface
  double precision,intent(in)::spacing_lat_interface,spacing_lon_interface
  double precision,dimension(nlat_interface,nlon_interface)::latlon_interface


!auxiliary variables
  integer ::ilat,ilon
  integer ::icornerlat,icornerlon
  double precision ::dtopo_dlat,dtopo_dlon
  double precision ::dlat,dlon
  double precision ::lon_adjusted

! make the longitude range consistent with the input latlon_topo file.
! Thus, converted into the range -180deg -180deg or 0-360deg.
  if(orig_lon_interface.lt.0.0.and.lon.gt.180.0) then
     lon_adjusted=-360.0+lon
  else
     lon_adjusted=lon
  end if
  

  icornerlat=int((lat-orig_lat_interface)/spacing_lat_interface)+1
  icornerlon=int((lon_adjusted-orig_lon_interface)/spacing_lon_interface)+1
!  print *,'icor',icornerlat,icornerlon,orig_lat_interface,orig_lon_interface
  if(icornerlat<1.or.icornerlat>nlat_interface-1.or.icornerlon<1.or. &
     icornerlon>nlon_interface-1) then
     print *,'icor',lat,lon_adjusted,orig_lat_interface,orig_lon_interface
     stop 'the lat or lon_adjusted is out of the range provided!'
  end if
  dtopo_dlat=0.0
  dtopo_dlon=0.0
  dlat=lat-(orig_lat_interface+icornerlat*spacing_lat_interface)
  dlon=lon_adjusted-(orig_lon_interface+icornerlon*spacing_lon_interface)
  if(icornerlat<nlat_interface)  then
   dtopo_dlat=(latlon_interface(icornerlat+1,icornerlon)- latlon_interface(icornerlat+1,icornerlon))/ &
               spacing_lat_interface
  end if
  if(icornerlon<nlat_interface)  then
   dtopo_dlon=(latlon_interface(icornerlat,icornerlon+1)- latlon_interface(icornerlat,icornerlon+1))/ &
               spacing_lon_interface
  end if

  topography=latlon_interface(icornerlat,icornerlon)+dlat*dtopo_dlat+dlon*dtopo_dlon
                                        
end subroutine search_topo
