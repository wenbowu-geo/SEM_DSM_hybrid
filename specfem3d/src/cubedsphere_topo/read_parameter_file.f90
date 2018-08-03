!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine read_parameter_file(LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK,CUBED_SPHERE_PROJECTION, &
        NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE,LOCAL_PATH,SUPPRESS_UTM_PROJECTION,&
        INTERFACES_FILE,CAVITY_FILE,NSUBREGIONS,&
        USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

  implicit none

  include "constants.h"

  integer NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, UTM_MAX
  double precision LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

  logical SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH
  !WENBO
  logical CUBED_SPHERE_PROJECTION
  logical CREATE_ABAQUS_FILES,CREATE_DX_FILES

  integer NDOUBLINGS
  integer, dimension(2) :: ner_doublings

  character(len=256) LOCAL_PATH
  character(len=50) INTERFACES_FILE,CAVITY_FILE
!WENBO
! rotation matrix from Euler angles
    double precision, dimension(NDIM,NDIM) :: rotation_matrix
    
    double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! local variables
  integer NEX_MAX

  double precision DEPTH_BLOCK_KM!,RECORD_LENGTH_IN_SECONDS,hdur,minval_hdur

!  character(len=256) dummystring
  integer ierr
  integer, external :: err_occurred

! subregions parameters
  integer NSUBREGIONS
  integer ix_beg_region,ix_end_region,iy_beg_region,iy_end_region
  integer iz_beg_region,iz_end_region,imaterial_number

! material properties

  integer i,ireg,imat,idoubl

! open parameter file
  open(unit=IIN,file=MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES)) &
       //'Mesh_Par_file',status='old',action='read')

  call read_value_double_precision(IIN,IGNORE_JUNK,LATITUDE_MIN, 'mesher.LATITUDE_MIN')
  if(err_occurred() /= 0) return
  call read_value_double_precision(IIN,IGNORE_JUNK,LATITUDE_MAX, 'mesher.LATITUDE_MAX')
  if(err_occurred() /= 0) return
  call read_value_double_precision(IIN,IGNORE_JUNK,LONGITUDE_MIN, 'mesher.LONGITUDE_MIN')
  if(err_occurred() /= 0) return
  call read_value_double_precision(IIN,IGNORE_JUNK,LONGITUDE_MAX, 'mesher.LONGITUDE_MAX')
  if(err_occurred() /= 0) return
  call read_value_double_precision(IIN,IGNORE_JUNK,DEPTH_BLOCK_KM, 'mesher.DEPTH_BLOCK_KM')
  if(err_occurred() /= 0) return
  call read_value_integer(IIN,IGNORE_JUNK,UTM_PROJECTION_ZONE, 'mesher.UTM_PROJECTION_ZONE')
  if(err_occurred() /= 0) return
  call read_value_logical(IIN,IGNORE_JUNK,SUPPRESS_UTM_PROJECTION, 'mesher.SUPPRESS_UTM_PROJECTION')
  if(err_occurred() /= 0) return
  call read_value_logical(IIN,IGNORE_JUNK,CUBED_SPHERE_PROJECTION,'mesher.CUBED_SPHERE_PROJECTION_PROJECTION')
  if(err_occurred() /= 0) return

  call read_value_string(IIN,IGNORE_JUNK,INTERFACES_FILE, 'mesher.INTERFACES_FILE')
  if(err_occurred() /= 0) return

  call read_value_string(IIN,IGNORE_JUNK,CAVITY_FILE, 'mesher.INTERFACES_FILE')
  if(err_occurred() /= 0) return

  call read_value_integer(IIN,IGNORE_JUNK,NEX_XI, 'mesher.NEX_XI')
  if(err_occurred() /= 0) return
  call read_value_integer(IIN,IGNORE_JUNK,NEX_ETA, 'mesher.NEX_ETA')
  if(err_occurred() /= 0) return
  call read_value_integer(IIN,IGNORE_JUNK,NPROC_XI, 'mesher.NPROC_XI')
  if(err_occurred() /= 0) return
  call read_value_integer(IIN,IGNORE_JUNK,NPROC_ETA, 'mesher.NPROC_ETA')
  if(err_occurred() /= 0) return


  Z_DEPTH_BLOCK = - dabs(DEPTH_BLOCK_KM) * 1000.d0

! check that parameters computed are consistent
if(.not.CUBED_SPHERE_PROJECTION) then
  if(UTM_X_MIN >= UTM_X_MAX) stop 'horizontal dimension of UTM block incorrect'
  if(UTM_Y_MIN >= UTM_Y_MAX) stop 'vertical dimension of UTM block incorrect'
end if
! set time step and radial distribution of elements
! right distribution is determined based upon maximum value of NEX
  NEX_MAX = max(NEX_XI,NEX_ETA)

!WENBO
  if(.not.CUBED_SPHERE_PROJECTION) then
     UTM_MAX = max(UTM_Y_MAX-UTM_Y_MIN, UTM_X_MAX-UTM_X_MIN)/1000.0 ! in KM
  end if
!WENBO
  call read_value_logical(IIN,IGNORE_JUNK,USE_REGULAR_MESH, 'mesher.USE_REGULAR_MESH')
  if(err_occurred() /= 0) return
  call read_value_integer(IIN,IGNORE_JUNK,NDOUBLINGS, 'mesher.NDOUBLINGS')
  if(err_occurred() /= 0) return
  call read_value_integer(IIN,IGNORE_JUNK,ner_doublings(1), 'mesher.NZ_DOUGLING_1')
  if(err_occurred() /= 0) return
  call read_value_integer(IIN,IGNORE_JUNK,ner_doublings(2), 'mesher.NZ_DOUGLING_2')
  if(err_occurred() /= 0) return

  if(ner_doublings(1) < ner_doublings(2) .and. NDOUBLINGS == 2) then
    idoubl = ner_doublings(1)
    ner_doublings(1) = ner_doublings(2)
    ner_doublings(2) = idoubl
  end if



  close(IIN)

  end subroutine read_parameter_file
