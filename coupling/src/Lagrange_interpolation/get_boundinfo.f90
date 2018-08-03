subroutine get_boundinfo()
  use coupling_SEM_DSM_par
  implicit none

  integer ::ipack_local,global_ipack_old,ios
  integer ::ipoint,ipoint_read
!other variables
  character(len=256) ::old_file,new_dir


!  write(new_dir,"('./BOUND_FFT/iproc',i4.4)") myrank
!  call system('mkdir -p '//adjustl(trim(new_dir)))

  ipoint=0
  do ipack_local=1,npack
     global_ipack_old=ipack_local+ipack_start-1
     write(old_file,"('./SEM_input/iproc',i4.4,'/bound_info_pack',i6.6)") &
                   in_iproc(global_ipack_old),local_pack_id(global_ipack_old)
     open(unit=53,file=trim(old_file),action='read',form='formatted',status="old",iostat=ios)
     do ipoint_read=1,npoint_pack(global_ipack_old)
       ipoint=ipoint+1
       read(53,*)id_depth(ipoint),is_elastic(ipoint),is_acoustic(ipoint)
!       if(id_depth(ipoint).eq.1) then
!          id_depth(ipoint)=1
!       else if(id_depth(ipoint).eq.2) then
!          id_depth(ipoint)=3
!       else if(id_depth(ipoint).lt.44) then
!          id_depth(ipoint)=(id_depth(ipoint)-2)*4+3
!       else 
!          id_depth(ipoint)=169
!       end if

!       read(53,*)c(:,ipoint)
       read(53,*)normal_vector(:,ipoint)
       read(53,*)jacobian(ipoint)
       read(53,*)coord_bound(:,ipoint)
!       if(id_depth(ipoint).eq.1) &
!          write(*,'(A10,2I4,5E13.5)'),"idep1_",local_pack_id(global_ipack_old),ipoint,&
!                coord_bound(2,ipoint),coord_bound(3,ipoint),normal_vector(:,ipoint)


     end do
     close(53)

  end do



end subroutine get_boundinfo
