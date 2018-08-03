subroutine compute_pml_rotation_matrix(myrank)

use constants, only:PI,CPML_Z_ONLY,CPML_X_ONLY,CPML_XZ_ONLY,NGLLX,&
        NGLLY,NGLLZ,TINYVAL,TINYVAL_SNGL,CUSTOM_REAL,NDIM
use specfem_par, only:ibool,xstore,ystore,zstore
use SEMtoTele_par, only:CENTER_LONGITUDE_IN_DEGREES,&
        CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,&
        ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES
use pml_par, only: NSPEC_CPML,CPML_regions,CPML_to_spec,pml_rotation_matrix
implicit none
integer, intent(in) ::myrank

real(kind=CUSTOM_REAL),parameter::radtodeg=180.0/PI
real(kind=CUSTOM_REAL),parameter::degtorad=PI/180.0
real(kind=CUSTOM_REAL),dimension(NDIM,NDIM) :: rotation_matrix0,&
        rotation_matrix_xi_up,rotation_matrix_xi_down,rotation_matrix_eta_up,&
        rotation_matrix_eta_down,rotation_matrix_back_cubedsph
real(kind=CUSTOM_REAL),dimension(NDIM,NDIM) ::pml_rotation_matrix_zbot,&
        pml_rotation_matrix_ximin,pml_rotation_matrix_ximax,pml_rotation_matrix_etamin,&
        pml_rotation_matrix_etamax
real(kind=CUSTOM_REAL) ::xi_min,xi_max,eta_min,eta_max
real(kind=CUSTOM_REAL) ::r,theta,phi,lat,lon,x,y,z
real(kind=CUSTOM_REAL) ::x_cubedsph,y_cubedsph,z_cubedsph

integer ::ispec_pml,ispec
!double precision, dimension(NDIM2D,NGNOD2D,NGLLY,NGLLZ):: dershape2D_x
!double precision, dimension(NDIM2D,NGNOD2D,NGLLX,NGLLZ):: dershape2D_y
!double precision, dimension(NDIM2D,NGNOD2D,NGLLX,NGLLY):: dershape2D_bottom
!double precision, dimension(NDIM2D,NGNOD2D,NGLLX,NGLLY):: dershape2D_top
!real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY)          :: jacobian2Dw_face
!real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY)     ::normal_face
!real(kind=CUSTOM_REAL)::  r
!real(kind=CUSTOM_REAL), dimension(NDIM):: coord_normalized,vector3_orthog
!integer ::iface

!From global to cubed sphere coordinages associated with lat0 and lon0
call euler_angles(rotation_matrix0,CENTER_LONGITUDE_IN_DEGREES,&
             CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)
rotation_matrix_back_cubedsph=transpose(rotation_matrix0)
xi_max=ANGULAR_WIDTH_XI_IN_DEGREES*degtorad/2.0
xi_min=-xi_max
eta_max=ANGULAR_WIDTH_ETA_IN_DEGREES*degtorad/2.0
eta_min=-eta_max

!xi goes up-> y_local goes up, due to y_sphere = x*rt
rotation_matrix_xi_up(:,:)=0.0
rotation_matrix_xi_up(1,1)=1.0
rotation_matrix_xi_up(2,2)=cos(xi_max)
rotation_matrix_xi_up(2,3)=sin(xi_max)
rotation_matrix_xi_up(3,2)=-sin(xi_max)
rotation_matrix_xi_up(3,3)=cos(xi_max)

rotation_matrix_xi_down(:,:)=0.0
rotation_matrix_xi_down(1,1)=1.0
rotation_matrix_xi_down(2,2)=cos(xi_min)
rotation_matrix_xi_down(2,3)=sin(xi_min)
rotation_matrix_xi_down(3,2)=-sin(xi_min)
rotation_matrix_xi_down(3,3)=cos(xi_min)

!eta goes to positive-> x_local goes to negative, due to x_sphere = -y*rt
rotation_matrix_eta_up(:,:)=0.0
rotation_matrix_eta_up(1,1)=cos(eta_max)
rotation_matrix_eta_up(1,3)=-sin(eta_max)
rotation_matrix_eta_up(2,2)=1.0
rotation_matrix_eta_up(3,1)=sin(eta_max)
rotation_matrix_eta_up(3,3)=cos(eta_max)

rotation_matrix_eta_down(:,:)=0.0
rotation_matrix_eta_down(1,1)=cos(eta_min)
rotation_matrix_eta_down(1,3)=-sin(eta_min)
rotation_matrix_eta_down(2,2)=1.0
rotation_matrix_eta_down(3,1)=sin(eta_min)
rotation_matrix_eta_down(3,3)=cos(eta_min)


pml_rotation_matrix_ximax=matmul(rotation_matrix0,rotation_matrix_xi_up)
pml_rotation_matrix_ximin=matmul(rotation_matrix0,rotation_matrix_xi_down)
pml_rotation_matrix_etamax=matmul(rotation_matrix0,rotation_matrix_eta_up)
pml_rotation_matrix_etamin=matmul(rotation_matrix0,rotation_matrix_eta_down)

if(myrank.eq.0) print *,"rot_matrix_ximin",pml_rotation_matrix_ximin(:,1),&
  pml_rotation_matrix_ximin(:,2),pml_rotation_matrix_ximin(:,3)
if(myrank.eq.0) print *,"rot_matrix_ximax",pml_rotation_matrix_ximax(:,1),&
        pml_rotation_matrix_ximax(:,2),pml_rotation_matrix_ximax(:,3)
if(myrank.eq.0) print *,"rot_matrix_etamin",pml_rotation_matrix_etamin(:,1),&
  pml_rotation_matrix_etamin(:,2),pml_rotation_matrix_etamin(:,3)
if(myrank.eq.0) print *,"rot_matrix_etamax",pml_rotation_matrix_etamax(:,1),&
  pml_rotation_matrix_etamax(:,2),pml_rotation_matrix_etamax(:,3)

do ispec_pml=1,NSPEC_CPML
  ispec = CPML_to_spec(ispec_pml)
    x=xstore(ibool((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec))
    y=ystore(ibool((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec))
    z=zstore(ibool((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec))

!to be fixed
  if(CPML_regions(ispec)==CPML_Z_ONLY) then
!xyz to lat and lon
    r=sqrt(x**2+y**2+z**2)
    theta=acos(z/r)
    phi=acos(x/(r*sin(theta))*(1-TINYVAL))
    if(y<TINYVAL_SNGL) then
       phi=-phi+2*PI
    end if
    lat=90.0-theta*radtodeg
    lon=phi*radtodeg
!cubed sphere coordinates associated with this element itself
     call euler_angles(pml_rotation_matrix_zbot,lon,lat,GAMMA_ROTATION_AZIMUTH)
     pml_rotation_matrix(:,:,ispec_pml)=pml_rotation_matrix_zbot(:,:)
     if(myrank.eq.0) print *,"rot_matrix_Z",pml_rotation_matrix(:,:,ispec_pml)
!to be fixed
  else 
     call xyz_backto_cubedsph(x,y,z,x_cubedsph,y_cubedsph,z_cubedsph,&
           ANGULAR_WIDTH_XI_IN_DEGREES, &
           ANGULAR_WIDTH_ETA_IN_DEGREES,rotation_matrix_back_cubedsph)
    if(CPML_regions(ispec)==CPML_X_ONLY.or.CPML_regions(ispec)==CPML_XZ_ONLY) then
       if(y_cubedsph.gt.0) then
        pml_rotation_matrix(:,:,ispec_pml)=pml_rotation_matrix_etamax(:,:)
       else
        pml_rotation_matrix(:,:,ispec_pml)=pml_rotation_matrix_etamin(:,:)
       end if
    else   ! Y XY YZ XYZ
       if(x_cubedsph.gt.0) then
        pml_rotation_matrix(:,:,ispec_pml)=pml_rotation_matrix_ximax(:,:)
       else
        pml_rotation_matrix(:,:,ispec_pml)=pml_rotation_matrix_ximin(:,:)
       end if
    end if
  end if


!  call get_jacobian_boundary_face(myrank,NSPEC_AB, &
!           xstore,ystore,zstore,ibool,NGLOB_AB,&
!           dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top,&
!           dble(wgllwgll_xy),dble(wgllwgll_xz),dble(wgllwgll_yz),&
!           ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)
!  if (NGNOD2D == 9.or.NGNOD2D == 5) then
!     x=xstore(ibool((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec))
!     y=ystore(ibool((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec))
!     z=zstore(ibool((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec))
!  else
!     x=xstore(ibool(1,1,1,ispec))
!     y=ystore(ibool(1,1,1,ispec))
!     z=zstore(ibool(1,1,1,ispec))
!  end if
!  r=sqrt(coord_normalized(1)**2+coord_normalized(2)**2+coord_normalized(3)**2)
!  coord_normalized(1)=x/r
!  coord_normalized(2)=y/r
!  coord_normalized(3)=z/r
!  vector3_orthog(1)=coord_normalized(2)*normal_face(3)-coord_normalized(3)*normal_face(2)
!  vector3_orthog(2)=coord_normalized(3)*normal_face(1)-coord_normalized(1)*normal_face(3)
!  vector3_orthog(3)=coord_normalized(1)*normal_face(3)-coord_normalized(2)*normal_face(1)
!  rotation_matrix(:,3,ispec_pml)=coord_normalized(:)
end do

end subroutine compute_pml_rotation_matrix
