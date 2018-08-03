!This subroutine is directly stolen from specfem3D/src/meshfem3D.

 subroutine read_value_integer(iunit,ignore_junk,value_to_read,name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  integer :: iunit
  integer :: value_to_read
  integer :: ier
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read

  call unused_string(name)
  ier = 0

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_integer

!--------------------

  subroutine read_value_dble_precision(iunit,ignore_junk,value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  integer :: iunit
  double precision :: value_to_read
  integer :: ier
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read

  call unused_string(name)
  ier = 0

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_dble_precision

!--------------------


  subroutine read_value_logical(iunit,ignore_junk,value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  logical :: value_to_read
  integer :: iunit
  integer :: ier
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read

  call unused_string(name)
  ier = 0

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_logical

!--------------------



  subroutine read_value_string(iunit,ignore_junk,value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  integer :: iunit
  character(len=*) :: value_to_read
  integer :: ier
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read

  call unused_string(name)
  ier = 0

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  value_to_read = string_read

  end subroutine read_value_string


  subroutine read_value_doubling_integer(iunit,ignore_junk,value_to_read,name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  integer :: iunit
  integer :: value_to_read
  integer :: ier
  character(len=*) :: name
  ! local parameters
  character(len=MAX_STRING_LEN) :: string_read
  integer :: index_equal_sign

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  ! checks if line contains name string
  if (index(string_read,trim(name)) > 0) then

    ! suppress leading junk (up to the first equal sign, included)
    index_equal_sign = index(string_read,'=')
    if (index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) &
      stop 'incorrect syntax detected in Par_file'

    string_read = string_read(index_equal_sign + 1:len_trim(string_read))

    ! suppress leading and trailing white spaces again, if any, after having
    ! suppressed the leading junk
    string_read = adjustl(string_read)
    string_read = string_read(1:len_trim(string_read))

    read(string_read,*,iostat=ier) value_to_read
  else
    ! returns an error
    ier = 1
    return
  endif

  end subroutine read_value_doubling_integer


!--------------------

  subroutine read_next_line(iunit,suppress_junk,string_read, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer :: iunit,ier
  logical :: suppress_junk
  character(len=MAX_STRING_LEN) :: string_read
  integer :: index_equal_sign

  ier = 0
  do
    read(unit=iunit,fmt="(a)",iostat=ier) string_read
    if (ier /= 0) stop 'Error while reading parameter file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text
! file coming from Windows/DOS)
    if (index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

! exit loop when we find the first line that is not a comment or a white line
    if (len_trim(string_read) == 0) cycle
    if (string_read(1:1) /= '#') exit

  enddo

! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

! suppress trailing comments, if any
  if (index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

  if (suppress_junk) then
! suppress leading junk (up to the first equal sign, included)
     index_equal_sign = index(string_read,'=')
     if (index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) &
      stop 'incorrect syntax detected in Par_file'

     string_read = string_read(index_equal_sign + 1:len_trim(string_read))
  endif

! suppress leading and trailing white spaces again, if any, after having
! suppressed the leading junk
  string_read = adjustl(string_read)
  string_read = string_read(1:len_trim(string_read))

  end subroutine read_next_line


!--------------------

! dummy subroutine to avoid warnings about variable not used in other
! subroutines
  subroutine unused_string(s)

  character(len=*) s

  if (len(s) == 1) continue

  end subroutine unused_string
