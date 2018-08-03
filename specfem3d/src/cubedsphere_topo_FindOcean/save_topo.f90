subroutine save_topo(nx_interface,ny_interface,xy_interface,interface_file,&
          orig_x_interface,orig_y_interface,&
          spacing_x_interface,spacing_y_interface)
  implicit none
  
  include "constants.h"

  double precision ::orig_x_interface,orig_y_interface
  double precision ::spacing_x_interface,spacing_y_interface
  integer ::nx_interface,ny_interface
  double precision::xy_interface(nx_interface,ny_interface)
  character(len=50)::interface_file

!auxiliary variables
  integer ::ix,iy

  open(unit=15,file=MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES))//interface_file, &
       status='unknown',action='write',form='formatted')
  write(15,*) nx_interface,ny_interface
  write(15,*) orig_x_interface,orig_y_interface
  write(15,*) spacing_x_interface,spacing_y_interface
  do iy=1,ny_interface
    do ix=1,nx_interface
        write(15,*) xy_interface(ix,iy)
    end do
  end do
  close(15)


end subroutine save_topo
