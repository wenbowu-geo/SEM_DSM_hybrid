subroutine  find_slab(x,y,z,x_cubedsph,y_cubedsph,z_cubedsph,vp,vs,rho)
use slabdepth_model
double precision:x,y,z,x_cubedsph,y_cubedsph,z_cubedsph
real(kind=CUSTOM_REAL), intent(out) ::vp,vs,rho

double precision::theta,phi,lat,lon,radi,depth
integer ::ilat_find,ilon_find,ilat,ilon
double precision::perturb
double precision::dist,dist_min
double precision, parameter::slab_thickness_toplayer=80000.0
double precision, parameter::slab_thickness_transition=20000.0
double precision, parameter::max_perterb=0.02
double precision, parameter::dlat_searchrange=4.0
double precision, parameter::dlon_searchrange=4.0

  radi=dsqrt(x**2+y**2+z**2)
  theta=dacos(z/radi)
  if(y<0.0) then
       phi=2*PI-dacos(x*(1.d0-1.e-10)/(radi*dsin(theta)))
  else
       phi=dacos(x*(1.d0-1.e-10)/(radi*dsin(theta)))
  end if
  lat= 180/PI*(PI/2-theta)
  lon=phi/PI*180.0
  depth=R_EARTH_SURF-radi

ilat_find=(lat-lat0_slabmodel)/dlat_slabmodel+1
ilon_find=(lon-lon0_slabmodel)/dlon_slabmodel+1

!check whether slab exists underground
if(slabdepth(ilon_find,ilat_find)<0.0) then
   perturb=0.0
else if(depth<slabdepth(ilon_find,ilat_find))
else 
!find the closest point on slab interface.
   dist_min=1.e10
   do ilat=1,nlat_slabmodel
     do ilon=1,nlon_slabmodel
        dist=(x-slabx(ilat,ilon))**2+(y-slaby(ilat,ilon))**2+(z-slabz(ilat,ilon))**2
        if(dist_min>dist) then
           dist_min=dist
           ilat_find=ilat
           ilon_find=ilon
        end if
     end do
   end do
   if(dist_min<slab_thickness_toplayer) then
      perturb=max_perterb*slabtaper(ilon,ilat)
   else if(dist_min<slab_thickness_toplayer+slab_thickness_transition) then
      perturb=max_perterb*(1-(dist_min-slab_thickness_toplayer)/slab_thickness_transition)*slabtaper(ilon,ilat)
   else
      perturb=0.0
   end if
end if

!add taper on the boundaries of simulation box
vp=vp*(1.0+perturb)
vs=vs*(1.0+perturb)
rho=rho*(1.0+perturb)

end subroutine
