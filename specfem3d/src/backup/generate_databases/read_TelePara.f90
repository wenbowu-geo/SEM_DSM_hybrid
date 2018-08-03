subroutine read_parameter_SEMtoTele(myrank)

!include "constants.h"

use constants
use SEMtoTele_par
implicit none


integer myrank
!integer::nx_TopoTaper,ny_TopoTaper,nx_notopo,ny_notopo
!integer::Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,Tele_irTop,Tele_irBot
!integer::NEX_XI,NEX_ETA,NER,NDOUBLINGS
!logical USE_REGULAR_MESH
!integer::ner_doublings(2)
!integer::myrank
!integer::Imain_SEMtoTele

integer :: ier


! open parameter file
 open(unit=IIN,file=MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES)) &
       //'SEMtoTele_Par_file',status='old',action='read')
 call read_value_double_precision_tele(IIN,IGNORE_JUNK,ANGULAR_WIDTH_XI_IN_DEGREES,'SEMtoTele.ANGULAR_WIDTH_XI_IN_DEGREES',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter ANGULAR_WIDTH_XI_IN_DEGREES'

 call read_value_double_precision_tele(IIN,IGNORE_JUNK,ANGULAR_WIDTH_ETA_IN_DEGREES,'SEMtoTele.ANGULAR_WIDTH_ETA_IN_DEGREES',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter ANGULAR_WIDTH_ETA_IN_DEGREES'

 call read_value_double_precision_tele(IIN,IGNORE_JUNK,CENTER_LATITUDE_IN_DEGREES,'SEMtoTele.CENTER_LATITUDE_IN_DEGREES',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter CENTER_LATITUDE_IN_DEGREES'

 call read_value_double_precision_tele(IIN,IGNORE_JUNK,CENTER_LONGITUDE_IN_DEGREES,'SEMtoTele.CENTER_LONGITUDE_IN_DEGREES',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter CENTER_LONGITUDE_IN_DEGREES'

 call read_value_double_precision_tele(IIN,IGNORE_JUNK,GAMMA_ROTATION_AZIMUTH,'SEMtoTele.GAMMA_ROTATION_AZIMUTH',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter GAMMA_ROTATION_AZIMUTH'

! call read_value_double_precision_tele(IIN,IGNORE_JUNK,DEPTH_BLOCK_KM_tmp,'SEMtoTele.DEPTH_BLOCK_KM',ier)
! if(ier /= 0) stop 'Error reading SEMtoTele paramete DEPTH_BLOCK_KM'

! Z_DEPTH_BLOCK = - dabs(DEPTH_BLOCK_KM_tmp) * 1000.d0


 call read_value_integer_tele(IIN,IGNORE_JUNK,NX_TOPOTAPER, 'SEMtoTele.nx_TopographyTaper',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter NX_TOPOTAPER'

 call read_value_integer_tele(IIN,IGNORE_JUNK,NY_TOPOTAPER, 'SEMtoTele.ny_TopographyTaper',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter NY_TOPOTAPER'

 call read_value_integer_tele(IIN,IGNORE_JUNK,NX_NOTOPO, 'SEMtoTele.nx_Notopograpgy',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter NX_NOTOPO'

 call read_value_integer_tele(IIN,IGNORE_JUNK,NY_NOTOPO, 'SEMtoTele.ny_Notopograpgy',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter NY_NOTOPO'

 call read_value_integer_tele(IIN,IGNORE_JUNK,TELE_IXLOW, 'SEMtoTele.ix_LowBound',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IXLOW'

 call read_value_integer_tele(IIN,IGNORE_JUNK,TELE_IXHIGH, 'SEMtoTele.ix_HighBound',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IXHIGH'

 call read_value_integer_tele(IIN,IGNORE_JUNK,TELE_IYLOW, 'SEMtoTele.iy_LowBound',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IYLOW'

 call read_value_integer_tele(IIN,IGNORE_JUNK,TELE_IYHIGH, 'SEMtoTele.iy_HighBound',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IYHIGH'

 call read_value_integer_tele(IIN,IGNORE_JUNK,TELE_IRTOP, 'SEMtoTele.ir_Top',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IRTOP'

 call read_value_integer_tele(IIN,IGNORE_JUNK,TELE_IRBOT, 'SEMtoTele.ir_Bound',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IRBOT'

 call read_value_logical_tele(IIN,IGNORE_JUNK,LOW_RESOLUTION,'SEMtoTele.LOW_RESOLUTION',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter LOW_RESOLUTION'

 call read_value_integer_tele(IIN,IGNORE_JUNK,NSTEP_BETWEEN_OUTPUTBOUND,'SEMtoTele.NSTEP_BETWEEN_OUTPUTBOUND',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter NSTEP_BETWEEN_OUTPUTBOUND'

 call read_value_integer_tele(IIN,IGNORE_JUNK,DECIMATE_COUPLING,'SEMtoTele.DECIMATE_COUPLING',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter DECIMATE_COUPLING'

 call read_value_integer_tele(IIN,IGNORE_JUNK,NPOINTS_PER_PACK,'SEMtoTele.NPOINTS_PER_PACK',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter NPOINTS_PER_PACK'

! close parameter file
 close(IIN)
 call euler_angles_back_cubedsph(rotation_matrix_back_cubedsph,CENTER_LONGITUDE_IN_DEGREES,&
                                 CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)



end subroutine read_parameter_SEMtoTele



! read values from parameter file, ignoring white lines and comments

  subroutine read_value_integer_tele(iunit,ignore_junk,value_to_read, name,ier)

  implicit none

  logical ignore_junk
  integer iunit
  integer value_to_read
  character(len=*) name
  character(len=100) string_read
  integer :: ier

  call unused_string_tele(name)

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  read(string_read,*) value_to_read

  end subroutine read_value_integer_tele

!--------------------

  subroutine read_value_double_precision_tele(iunit,ignore_junk,value_to_read, name,ier)

  implicit none

  logical ignore_junk
  integer iunit
  double precision value_to_read
  character(len=*) name
  character(len=100) string_read
  integer :: ier

  call unused_string_tele(name)

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  read(string_read,*) value_to_read

  end subroutine read_value_double_precision_tele

!--------------------

  subroutine read_value_logical_tele(iunit,ignore_junk,value_to_read, name,ier)

  implicit none

  logical ignore_junk
  logical value_to_read
  integer iunit
  character(len=*) name
  character(len=100) string_read
  integer :: ier

  call unused_string_tele(name)

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  read(string_read,*) value_to_read

  end subroutine read_value_logical_tele

!--------------------

  subroutine read_value_string_tele(iunit,ignore_junk,value_to_read, name,ier)

  implicit none

  logical ignore_junk
  integer iunit
  character(len=*) value_to_read
  character(len=*) name
  character(len=100) string_read
  integer :: ier

  call unused_string_tele(name)

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  value_to_read = string_read

  end subroutine read_value_string_tele

!--------------------

!--------------------

  subroutine read_next_line(iunit,suppress_junk,string_read,ier)

  implicit none

  include "constants.h"


  logical suppress_junk
  character(len=100) string_read
  integer index_equal_sign,ier,iunit

  ier = 0
  do
    read(unit=iunit,fmt="(a100)",iostat=ier) string_read
    if(ier /= 0) stop 'error while reading parameter file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text
! file coming from Windows/DOS)
    if(index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

! exit loop when we find the first line that is not a comment or a white line
    if(len_trim(string_read) == 0) cycle
    if(string_read(1:1) /= '#') exit

  enddo

! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

! suppress trailing comments, if any
  if(index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

  if(suppress_junk) then
! suppress leading junk (up to the first equal sign, included)
     index_equal_sign = index(string_read,'=')
     if(index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) stop 'incorrect syntax detected in Mesh_Par_file'
     string_read = string_read(index_equal_sign + 1:len_trim(string_read))
  end if

! suppress leading and trailing white spaces again, if any, after having
! suppressed the leading junk
  string_read = adjustl(string_read)
  string_read = string_read(1:len_trim(string_read))

  end subroutine read_next_line

!--------------------

!--------------------

  integer function err_occurred_tele()

  err_occurred_tele = 0

  end function err_occurred_tele

!--------------------

! dummy subroutine to avoid warnings about variable not used in other
! subroutines
  subroutine unused_string_tele(s)

  character(len=*) s

  if (len(s) == 1) continue

  end subroutine unused_string_tele
