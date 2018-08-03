subroutine save_depth_id_innerbound(myrank,NPROC,Nele_Bound,&
                     ispec_is_elastic,&
                     ispec_is_acoustic,nspec,ibool,nglob,xstore_dummy,&
                     ystore_dummy,zstore_dummy,media_type,wzgll,&
                     IMAIN_OUTPUT,&
                     prname,LOCAL_PATH)
  use tomography, only:rmin_structure_zone,n_structure_zone
  use SEMtoTele_par, only:Bound_Info,LOW_RESOLUTION,id_depth_bound,id_dist_bound
  use constants
  implicit none

!  include "constants.h"

  integer ::myrank,NPROC,media_type
  integer :: Nele_Bound
  integer :: IMAIN_OUTPUT

  integer nspec,nglob
! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy
  double precision, dimension(NGLLZ) ::wzgll

  logical, dimension(nspec) :: ispec_is_acoustic,ispec_is_elastic
  character(len=256) prname,LOCAL_PATH


!other parameters
  integer ::npoint_bound,ipoint,ipoint_iface,ipoint_temp,ipoint_allmedia
  integer ::max_ndep_global,ndepth_table,&
            idepth,ndist_table,idist
  integer ::num_faces_this_mediatype
  integer ::size_depth_bound,size_dist_bound,max_ndist_iproc,max_ndist_allproc
  integer ::max_ndist_global
  integer, allocatable,dimension(:) ::ndist
  real(kind=CUSTOM_REAL) ::r_thisdepth,r_within_mesh
  integer ::izone,idiscontinuity,iglob_within_mesh
  integer,allocatable,dimension(:) ::depth_id,dist_id
  integer,allocatable,dimension(:) ::iface_this_mediatype
  real(kind=CUSTOM_REAL) ::dist_rad
  real(kind=CUSTOM_REAL), allocatable, dimension(:) ::depth_bound,dist_bound
  real(kind=CUSTOM_REAL), allocatable, dimension(:) ::depth_table_global,dist_table_global
  integer, allocatable, dimension(:) ::izone_dep_table,move_idepth,discon_idepth
  integer ::ndiscon_found_deptable
  integer ::ndep_discon_acount
  real(kind=CUSTOM_REAL), allocatable, dimension(:) ::dep_table_discon_acount
  integer, allocatable, dimension(:) ::izone_dep_discon_acount
  integer,dimension(:),allocatable ::idepth_is_discont
  real(kind=CUSTOM_REAL), allocatable, dimension(:) ::dist_thisdepth
  integer,allocatable,dimension(:) ::dist_id_thisdepth
  double precision ::x,y,z,r
  
  double precision,dimension(3) ::coord_sta,coord_bound
  integer ::i,j,k,igll,iface,iface_temp,ispec,iglob
  integer ::imax,imin,jmax,jmin,kmax,kmin,di,dj,dk

!parameters for sending and receiving
  integer::iproc,req_send,req_recv

!Tele-station
  double precision :: lat_sta(1),long_sta(1),depth_sta(1)
  double precision ::t_cmt(1),hdur(1)
  double precision ::theta_sta,phi_sta
  double precision moment_tensor(6,1)

  double precision sec


  character(len=256) file_name,clean_LOCAL_PATH
  character(len=21) :: elastic_file
  character(len=22) :: acoustic_file
  data elastic_file /'id_depth_elas_innerbound.info'/
  data acoustic_file /'id_depth_acous_innerbound.info'/


  num_faces_this_mediatype=0
  do iface=1,Nele_Bound
     ispec = Bound_Info(iface)%ispec_bound
!    if(DEBUG_COUPLING) print *,'iface',iface,ispec_is_elastic(ispec)
     if(ispec_is_acoustic(ispec).and.media_type.eq.2) then
        num_faces_this_mediatype=num_faces_this_mediatype+1
     else if(ispec_is_elastic(ispec).and.media_type.eq.1) then
        num_faces_this_mediatype=num_faces_this_mediatype+1
     end if
  end do
  if(DEBUG_COUPLING) print *,'num_faces_this_mediatype',num_faces_this_mediatype
  if(num_faces_this_mediatype.gt.0) then
     allocate(iface_this_mediatype(num_faces_this_mediatype))
  else 
       if(media_type.eq.1) then
         file_name = prname(1:len_trim(prname))//elastic_file(1:len_trim(elastic_file))
       else
         file_name = prname(1:len_trim(prname))//acoustic_file(1:len_trim(acoustic_file))
       end if

     open(unit=16,file=file_name(1:len_trim(file_name)),& 
          status='unknown',action='write',form='formatted')

     write(16,*) 0

     close(16)
  end if
  if(DEBUG_COUPLING) print *,'num_faces_this_mediatype',num_faces_this_mediatype,media_type
  iface_temp=0
  do iface=1,Nele_Bound
     ispec = Bound_Info(iface)%ispec_bound
     if(ispec_is_acoustic(ispec).and.media_type.eq.2) then
        iface_temp=iface_temp+1
        iface_this_mediatype(iface_temp)=iface
        if(DEBUG_COUPLING) print *,'acoustic ispec',iface
     else if(ispec_is_elastic(ispec).and.media_type.eq.1) then
        iface_temp=iface_temp+1
        iface_this_mediatype(iface_temp)=iface
     else if(.not.ispec_is_acoustic(ispec).and..not.ispec_is_elastic(ispec)) then
        print *,'media type',ispec_is_acoustic(ispec),ispec_is_elastic(ispec),media_type,&
               .not.ispec_is_acoustic(ispec).and..not.ispec_is_acoustic(ispec)
        call exit_MPI(myrank,'unknown or wrong media type')
     end if
  end do


  if(LOW_RESOLUTION) then
      npoint_bound=num_faces_this_mediatype
  else
      npoint_bound=NGLLX*NGLLY*num_faces_this_mediatype
  end if
  if(DEBUG_COUPLING) print *,'depth_list begin',myrank,npoint_bound

  if(npoint_bound.gt.0) then
     size_depth_bound=npoint_bound
     size_dist_bound=npoint_bound
  else
     size_depth_bound=1
     size_dist_bound=1
  end if
  allocate(depth_bound(size_depth_bound))
  allocate(depth_id(size_depth_bound))
  allocate(dist_bound(size_dist_bound))
  allocate(dist_id(size_dist_bound))


!  if(myrank.eq.0) then
!   allocate(ndepth_eachproc(NPROC))
!   ndepth_eachproc(:)=0
!  end if
  
  if(DEBUG_COUPLING) print *,'allocate finished',myrank


  lat_sta(1)=0.0
  long_sta(1)=152.75 
  theta_sta = PI/2.0d0 - lat_sta(1)*PI/180.0d0
  phi_sta = long_sta(1)*PI/180.0d0
  coord_sta(1)=dsin(theta_sta)*dcos(phi_sta)
  coord_sta(2)=dsin(theta_sta)*dsin(phi_sta)
  coord_sta(3)=dcos(theta_sta)
  if(DEBUG_COUPLING) print *,'sta done'

  ipoint=0
  do iface=1,num_faces_this_mediatype
     ispec = Bound_Info(iface_this_mediatype(iface))%ispec_bound
     imin=Bound_Info(iface)%istart;imax=Bound_Info(iface)%iend
     jmin=Bound_Info(iface)%jstart;jmax=Bound_Info(iface)%jend
     kmin=Bound_Info(iface)%kstart;kmax=Bound_Info(iface)%kend
     di=1;dj=1;dk=1

     ipoint_iface=0
     do k=kmin,kmax,dk
      do j=jmin,jmax,dj
       do i=imin,imax,di
           if(LOW_RESOLUTION) then
               ipoint=iface
           else
               ipoint_iface=ipoint_iface+1
               ipoint=(iface-1)*NGLLX*NGLLY+ipoint_iface
           end if


           iglob=ibool(i,j,k,ispec)
           x=xstore_dummy(iglob)
           y=ystore_dummy(iglob)
           z=zstore_dummy(iglob)
           r=dsqrt(x**2+y**2+z**2)
           depth_bound(ipoint)=R_EARTH_SURF-r        
           coord_bound(1)=x/r;coord_bound(2)=y/r;coord_bound(3)=z/r

           dist_rad=dacos(dot_product(coord_bound,coord_sta)*(1.d0-TINYVAL))
!WENBO
!       dist_bound(ipoint) = r*dist_rad
           dist_bound(ipoint) = dist_rad*180.d0/PI

           if(depth_bound(ipoint)<100.and.DEBUG_COUPLING) &
              print *,'check_dist_reslult',iface,igll,x,y,z,dist_rad*180.0/3.1415926,myrank
       end do !i
      end do !j
     end do !k
  end do

if(DEBUG_COUPLING) print *,'MAX_NDEPTH_TABLE',MAX_NDEPTH_TABLE
max_ndep_global=MAX_NDEPTH_TABLE
allocate(depth_table_global(max_ndep_global))

if(DEBUG_COUPLING) print *,'find_global_id'
!construct depth table and get depth id
call find_global_id(depth_bound,size_depth_bound,npoint_bound,depth_tolerence,depth_id,&
                      depth_table_global,ndepth_table,max_ndep_global,&
                      NPROC,myrank)

if(myrank.eq.0.and.DEBUG_COUPLING) print *,'ndepth_table',ndepth_table,depth_table_global(:)

ndiscon_found_deptable=0
if(ndepth_table.gt.0) then


 if(media_type.eq.1) then
     file_name = trim(LOCAL_PATH(1:len_trim(LOCAL_PATH)))//'/dist_table_elastic_innerbound'
 else
     file_name = trim(LOCAL_PATH(1:len_trim(LOCAL_PATH)))//'/dist_table_acoustic_innerbound'
 end if
 if(myrank.eq.0) &
     open(113,file=file_name(1:len_trim(file_name)),status='unknown',action='write',form='formatted')


  allocate(ndist(ndepth_table))
  ndist(:)=0
  dist_id(:)=0
  max_ndist_iproc=1
  ipoint_temp=0
  do idepth=1,ndepth_table
    do ipoint=1,npoint_bound
        if(depth_id(ipoint).eq.idepth) then
             ndist(idepth)=ndist(idepth)+1
             ipoint_temp=ipoint_temp+1
        end if
    end do
    if(ndist(idepth).gt.max_ndist_iproc) max_ndist_iproc=ndist(idepth)
  end do
  if(npoint_bound.gt.0.and.ipoint_temp.ne.npoint_bound) then
        print *,myrank,"Error counting coupling points",ipoint_temp,npoint_bound
        call exit_MPI(myrank,'Error in counting coupling points')
  end if

  allocate(dist_thisdepth(max_ndist_iproc))
  allocate(dist_id_thisdepth(max_ndist_iproc))
  call  max_all_i(max_ndist_iproc,max_ndist_allproc)

  if(myrank.eq.0) then
    max_ndist_global=NPROC*max_ndist_allproc
  else 
    max_ndist_global=1
  end if
  allocate(dist_table_global(max_ndist_global))


  if(DEBUG_COUPLING) print *,'discon_start',myrank


!take the conscontinuity into account
  allocate(move_idepth(ndepth_table))
  allocate(idepth_is_discont(ndepth_table))
 if(myrank.eq.0) then
  allocate(izone_dep_table(ndepth_table))
  allocate(discon_idepth(ndepth_table))
  discon_idepth(:)=0
  izone_dep_table(:)=0
  move_idepth(:)=0
  ndiscon_found_deptable=0
  idepth_is_discont(:)=0
  do idepth=1,ndepth_table
    r_thisdepth=R_EARTH_SURF-depth_table_global(idepth)
    if(idepth.eq.1) r_thisdepth=r_thisdepth-0.5*depth_tolerence
    if(idepth.eq.ndepth_table) r_thisdepth=r_thisdepth+0.5*depth_tolerence
    do izone=1,n_structure_zone
      if(rmin_structure_zone(izone).le.r_thisdepth) then
         izone_dep_table(idepth)=izone
      end if
      if(idepth.eq.25.and.izone.eq.10.and.DEBUG_COUPLING) print *,'discon detect',&
            rmin_structure_zone(izone),r_thisdepth,depth_table_global(:)
      if(dabs(rmin_structure_zone(izone)-r_thisdepth).lt.0.5*depth_tolerence.and.&
         idepth.gt.1.and.idepth.lt.ndepth_table) then
         idepth_is_discont(idepth)=1
         ndiscon_found_deptable=ndiscon_found_deptable+1
         discon_idepth(ndiscon_found_deptable)=idepth
         izone_dep_table(idepth)=izone
         move_idepth(idepth:ndepth_table)=move_idepth(idepth:ndepth_table)+1
      end if
    end do
  end do
  if(izone_dep_table(1)-izone_dep_table(ndepth_table).ne.ndiscon_found_deptable.and.&
     .not.LOW_RESOLUTION) then
      print *,'check_izone',izone_dep_table(1),izone_dep_table(ndepth_table),&
            ndepth_table,ndiscon_found_deptable
      stop 'Error of discontinuity setting or depth table (innerbound)'
  end if
  ndep_discon_acount=ndepth_table+ndiscon_found_deptable

  allocate(dep_table_discon_acount(ndep_discon_acount))
  allocate(izone_dep_discon_acount(ndep_discon_acount))
  do idepth=1,ndepth_table
    dep_table_discon_acount(idepth+move_idepth(idepth))=depth_table_global(idepth)
    izone_dep_discon_acount(idepth+move_idepth(idepth))=izone_dep_table(idepth)
  end do


  do idiscontinuity=1,ndiscon_found_deptable
     dep_table_discon_acount(discon_idepth(idiscontinuity)+idiscontinuity-1)= &
          depth_table_global(discon_idepth(idiscontinuity))
     izone_dep_discon_acount(discon_idepth(idiscontinuity)+idiscontinuity-1)= &
          izone_dep_discon_acount(discon_idepth(idiscontinuity))
     izone_dep_discon_acount(discon_idepth(idiscontinuity)+idiscontinuity)= &
          izone_dep_discon_acount(discon_idepth(idiscontinuity))+1
  end do
 end if



  if(myrank.eq.0) then
!send ndep_discon_acount, discon_idepth and  arrays
    do iproc=1,NPROC-1

           req_send=iproc+1000
           call isend_i(idepth_is_discont(1),ndepth_table,&
                         iproc,1,req_send)
           call wait_req(req_send)
           req_send=iproc+1001
           call isend_i(move_idepth(1),ndepth_table,&
                         iproc,1,req_send)
           call wait_req(req_send)
    end do
  else
!receive the arrays
            req_recv=myrank+1000
            call irecv_i(idepth_is_discont(1),ndepth_table,&
                                        0,1,req_recv)
            call wait_req(req_recv)

            req_recv=myrank+1001
            call irecv_i(move_idepth(1),ndepth_table,&
                                        0,1,req_recv)
            call wait_req(req_recv)
  end if


 !construct distance table and get distance id
 do idepth=1,ndepth_table 
    if(DEBUG_COUPLING) print *,'construct distance table',idepth
    idist=0
    do ipoint=1,npoint_bound
      if(depth_id(ipoint).eq.idepth) then
          idist=idist+1
          dist_thisdepth(idist)=dist_bound(ipoint)
      end if
    end do

    call find_global_id(dist_thisdepth,max_ndist_iproc,ndist(idepth),dist_tolerence,&
                      dist_id_thisdepth,&
                      dist_table_global,ndist_table,max_ndist_global,&
                      NPROC,myrank)
!copy id
    idist=0
    do ipoint=1,npoint_bound
      if(depth_id(ipoint).eq.idepth) then
          idist=idist+1
          dist_id(ipoint)=dist_id_thisdepth(idist)
!          if(dist_id(ipoint).eq.0.and.DEBUG_COUPLING) print
!          *,"Error",dist_bound(ipoint),idist,ndist(idepth)
      end if
    end do
    if(idist.ne.ndist(idepth)) then
      print *,"Error in counting distnace",idist,ndist(idepth)
      call exit_MPI(myrank,'Error in counting distnace of coupling points.')
    end if


!save table
  if(myrank.eq.0) then
    r=R_EARTH_SURF-depth_table_global(idepth)
    write(113,*)ndist_table
    do idist=1,ndist_table
!      write(113,*) dist_table_global(idist)/r*(180.d0/PI)
      write(113,*) dist_table_global(idist)
    end do
    if(idepth_is_discont(idepth).eq.1) then
      write(113,*)ndist_table
      do idist=1,ndist_table
         write(113,*) dist_table_global(idist)
      end do

    end if

  end if
 end do  ! do idepth
 if(myrank.eq.0) close(113)

  if(DEBUG_COUPLING) print *,'move_idepth',move_idepth(:),idepth_is_discont(:)
 

 if(.not.LOW_RESOLUTION) then
  ipoint=0
  do iface=1,num_faces_this_mediatype
     ispec = Bound_Info(iface_this_mediatype(iface))%ispec_bound
     imin=Bound_Info(iface)%istart;imax=Bound_Info(iface)%iend
     jmin=Bound_Info(iface)%jstart;jmax=Bound_Info(iface)%jend
     kmin=Bound_Info(iface)%kstart;kmax=Bound_Info(iface)%kend
     di=1;dj=1;dk=1

     ipoint_iface=0
     do k=kmin,kmax,dk
      do j=jmin,jmax,dj
       do i=imin,imax,di
           if(LOW_RESOLUTION) then
               ipoint=iface
           else
               ipoint_iface=ipoint_iface+1
               ipoint=(iface-1)*NGLLX*NGLLY+ipoint_iface
           end if
          if(idepth_is_discont(depth_id(ipoint)).eq.1) then
            ispec=Bound_Info(iface_this_mediatype(iface))%ispec_bound 
            iglob=ibool(i,j,k,ispec)
            x=xstore_dummy(iglob)
            y=ystore_dummy(iglob)
            z=zstore_dummy(iglob)
            r=dsqrt(x**2+y**2+z**2)

            iglob_within_mesh=ibool(2,2,2,ispec)
            x=xstore_dummy(iglob_within_mesh)
            y=ystore_dummy(iglob_within_mesh)
            z=zstore_dummy(iglob_within_mesh)
            r_within_mesh=dsqrt(x**2+y**2+z**2)


            if(i.lt.NGLLX.and.i.gt.1 .and. j.lt.NGLLY.and.j.gt.1 &
              .and.k.lt.NGLLZ.and.k.gt.1 ) &
              stop 'Error, discontinuity is within this mesh!'

            if(r_within_mesh.lt.r) then
              depth_id(ipoint)=depth_id(ipoint)+move_idepth(depth_id(ipoint))
            else
              depth_id(ipoint)=depth_id(ipoint)+move_idepth(depth_id(ipoint)-1)
            end if
          else
              depth_id(ipoint)=depth_id(ipoint)+move_idepth(depth_id(ipoint))
          end if
       end do !do i
      end do ! do j
     end do !do k
   end do  !iface
  end if !LOW_RESOLUTION


  if(media_type.eq.1) then
     file_name = prname(1:len_trim(prname))//elastic_file(1:len_trim(elastic_file))
  else
     file_name = prname(1:len_trim(prname))//acoustic_file(1:len_trim(acoustic_file))
  end if


  open(unit=117,file=file_name(1:len_trim(file_name)),status='unknown',action='write',form='formatted')

  write(117,*) num_faces_this_mediatype
  if(DEBUG_COUPLING) print *,'num_faces',myrank,num_faces_this_mediatype,media_type
  ipoint=0
  do iface=1,num_faces_this_mediatype
     ispec = Bound_Info(iface_this_mediatype(iface))%ispec_bound
     imin=Bound_Info(iface)%istart;imax=Bound_Info(iface)%iend
     jmin=Bound_Info(iface)%jstart;jmax=Bound_Info(iface)%jend
     kmin=Bound_Info(iface)%kstart;kmax=Bound_Info(iface)%kend
     di=1;dj=1;dk=1

     ipoint_iface=0
     do k=kmin,kmax,dk
      do j=jmin,jmax,dj
       do i=imin,imax,di
          if(LOW_RESOLUTION) then
               ipoint=iface
               ipoint_allmedia=iface_this_mediatype(iface)
          else
               ipoint_iface=ipoint_iface+1
               ipoint=(iface-1)*NGLLX*NGLLY+ipoint_iface
               ipoint_allmedia=(iface_this_mediatype(iface)-1)*NGLLX*NGLLY+ipoint_iface
          end if
          write(117,*) iface_this_mediatype(iface),depth_id(ipoint),dist_id(ipoint),dist_bound(ipoint)
          id_depth_bound(ipoint_allmedia)=depth_id(ipoint)
          id_dist_bound(ipoint_allmedia)=dist_id(ipoint)
       end do
      end do
     end do
  end do

  close(117)


  deallocate(ndist)
  deallocate(dist_thisdepth)
  deallocate(dist_id_thisdepth)
  deallocate(dist_table_global)

end if

 !save depth table
  if(myrank.eq.0) then
     if(media_type.eq.1) then
       file_name = trim(LOCAL_PATH(1:len_trim(LOCAL_PATH)))//'/depth_table_elastic_innerbound'
     else
       file_name = trim(LOCAL_PATH(1:len_trim(LOCAL_PATH)))//'/depth_table_acoustic_innerbound'
     end if
     open(17,file=file_name(1:len_trim(file_name)),status='unknown',action='write',form='formatted')
!     write(17,*) NGLLZ
!     write(17,*) wzgll(1:NGLLZ)/2.d0
     write(17,*) 0.02,ndepth_table+ndiscon_found_deptable
     write(17,*) 0.005
     write(17,*) NGLLZ

     do idepth=1,ndepth_table
      write(17,*)depth_table_global(idepth)/1000.0,izone_dep_table(idepth)
      if(idepth_is_discont(idepth).eq.1) then
        write(17,*)depth_table_global(idepth)/1000.0,izone_dep_table(idepth)-1
      end if
     end do
     close(17)
  end if

 
!save depth and distance id
  


  if(num_faces_this_mediatype.gt.0) deallocate(iface_this_mediatype)

  deallocate(depth_bound)
  deallocate(depth_id)
  deallocate(depth_table_global)
  deallocate(dist_bound)
  deallocate(dist_id)
  if(ndepth_table.gt.0) then
    deallocate(move_idepth)
    deallocate(idepth_is_discont)
  end if

  if(myrank.eq.0.and.ndepth_table.gt.0) then
    deallocate(discon_idepth)
    deallocate(izone_dep_table)
    deallocate(izone_dep_discon_acount)
    deallocate(dep_table_discon_acount)
  end if



end subroutine



