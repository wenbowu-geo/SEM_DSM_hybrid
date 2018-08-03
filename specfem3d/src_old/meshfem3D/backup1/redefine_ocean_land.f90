subroutine redefine_ocean_land(xelm,yelm,zelm,imaterial_number,&
                            isubregion,nmeshregions,ir,ir1,ir2,iy,iy1,iy2,ix,ix1,ix2,&
                            iproc_xi,NPROC_XI,iproc_eta,NPROC_ETA,&
                            nx_notopo,ny_notopo)
  use constants
  use meshfem3D_par, only:myrank
  implicit none
!  include "constants.h"

 double precision xelm(NGNOD_EIGHT_CORNERS),yelm(NGNOD_EIGHT_CORNERS),zelm(NGNOD_EIGHT_CORNERS)
 integer imaterial_number
 integer isubregion,nmeshregions,ir,ir1,ir2,iy,iy1,iy2,ix,ix1,ix2
 integer iproc_xi,NPROC_XI,iproc_eta,NPROC_ETA
 integer nx_notopo,ny_notopo
 
 !local parameters
 double precision z_thismesh,y_thismesh,x_thismesh
 integer,save:: npx_interface_topo,npy_interface_topo
 double precision,save:: orig_x_interface_topo,orig_y_interface_topo
 double precision,save:: spacing_x_interface_topo,spacing_y_interface_topo
 double precision,dimension(:,:), allocatable,save:: interface_topo
 integer ix_tmp,iy_tmp,ix_find,iy_find
 integer ignod


 if(isubregion.eq.1.and.ir.eq.ir1.and.ix.eq.ix1.and.iy.eq.iy1) then

   open(unit=45,file=MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES))//'real_bathymetry_topography',status='old')
! call read_value_integer(45,DONT_IGNORE_JUNK,npx_interface_topo,'NINTERFACES') 
! call read_value_integer(45,DONT_IGNORE_JUNK,npy_interface_topo,'NINTERFACES') 
! call read_value_double_precision(45,DONT_IGNORE_JUNK,orig_x_interface_topo,'orig_x_interface_topo')
! call read_value_double_precision(45,DONT_IGNORE_JUNK,orig_y_interface_topo,'orig_y_interface_topo')
! call read_value_double_precision(45,DONT_IGNORE_JUNK,spacing_x_interface_topo,'spacing_x_interface_topo')
! call read_value_double_precision(45,DONT_IGNORE_JUNK,spacing_y_interface_topo,'spacing_y_interface_topo')
   read(45,*) npx_interface_topo,npy_interface_topo
   read(45,*) orig_x_interface_topo,orig_y_interface_topo
   read(45,*) spacing_x_interface_topo,spacing_y_interface_topo
   allocate(interface_topo(npx_interface_topo,npy_interface_topo))
 
   do iy_tmp=1,npy_interface_topo
     do ix_tmp=1,npx_interface_topo
!       call read_value_double_precision(45,DONT_IGNORE_JUNK,interface_topo(ix,iy),'Z_INTERFACE_TOPO')
         read(45,*) interface_topo(ix_tmp,iy_tmp)
     enddo
   enddo
   close(45)
 end if

 x_thismesh=0.d0
 y_thismesh=0.d0
 z_thismesh=0.d0
 do ignod=1,NGNOD_EIGHT_CORNERS
    x_thismesh=x_thismesh+xelm(ignod)
    y_thismesh=y_thismesh+yelm(ignod)
    z_thismesh=z_thismesh+zelm(ignod)
 end do
 x_thismesh=x_thismesh/NGNOD_EIGHT_CORNERS
 y_thismesh=y_thismesh/NGNOD_EIGHT_CORNERS
 z_thismesh=z_thismesh/NGNOD_EIGHT_CORNERS

 ix_find=int((x_thismesh - orig_x_interface_topo) / spacing_x_interface_topo) + 1
 iy_find = int((y_thismesh - orig_y_interface_topo) / spacing_y_interface_topo) + 1
! if(ix_find.gt.npx_interface_topo.or.ix_find.lt.1.or.iy_find.gt.npy_interface_topo.or.iy_find.lt.1) then
!    print *,'Error_detected',ix_find,iy_find,x_thismesh,y_thismesh
! end if

!box edges with no topograpny
 if((iproc_xi.eq.0.and.ix.le.nx_notopo).or.&
    (iproc_xi.eq.NPROC_XI-1.and.ix.ge.ix2-2*nx_notopo).or.&
    (iproc_eta.eq.0.and.iy.le.ny_notopo).or.&
    (iproc_eta.eq.NPROC_ETA-1.and.iy.ge.iy2-2*ny_notopo)) then 

!acoustic or elastic material depends on topography
!For some extreme cases where steep topography happens, the linearly interpolated
!z_thismesh in land might be higher than interfac_topo. Thus the element would
!be incorrectly classified to ocean. Using addtional requirement
!z_thismesh.lt.100 to exclude those cases.
   if(z_thismesh.gt.interface_topo(ix_find,iy_find).and.z_thismesh.lt.100) then
    !acoustic
    imaterial_number=1
   else
    !elastic
    imaterial_number=2
   end if

!pose elastic walls at four edges of box
   !elastic
!   imaterial_number=2
!For setting the solid walls at four edges of SEM box
!use imaterial_number=3 to mark the wall elements.
!   imaterial_number=3


!with topography
 else
   if(z_thismesh.gt.interface_topo(ix_find,iy_find).and.z_thismesh.lt.100) then
    !acoustic
    imaterial_number=1
   else
    !elastic
    imaterial_number=2
   end if
 end if
 
!debug
!if(myrank.eq.78.and.imaterial_number.eq.1) then
!   print *,'Error, water detected here,  it is wrong!'
!   print *,'location ix iy evelv',ix_find,iy_find,interface_topo(ix_find,iy_find),&
!        z_thismesh,zelm(:)
!end if

!WENBO
!to be deleted
!no ocean
!imaterial_number=2

 if(isubregion.eq.nmeshregions.and.ir.eq.ir2.and.ix.eq.ix2.and.iy.eq.iy2) then
   deallocate(interface_topo)
 end if
end subroutine redefine_ocean_land
