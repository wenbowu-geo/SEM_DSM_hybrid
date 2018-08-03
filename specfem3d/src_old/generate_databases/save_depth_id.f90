subroutine save_depth_id(myrank,NPROC,num_abs_boundary_faces,&
                     abs_boundary_ispec,abs_boundary_ijk,ispec_is_elastic,&
                     ispec_is_acoustic,nspec,ibool,nglob,xstore_dummy,&
                     ystore_dummy,zstore_dummy,media_type,wzgll,&
                     IMAIN_OUTPUT,prname,LOCAL_PATH)
  use tomography, only:rmin_structure_zone,n_structure_zone
  use constants
  implicit none

!  include "constants.h"

  integer ::myrank,NPROC,media_type
  integer :: num_abs_boundary_faces
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)
  integer :: IMAIN_OUTPUT

  integer nspec,nglob
! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy
  double precision, dimension(NGLLZ) ::wzgll

  logical, dimension(nspec) :: ispec_is_acoustic,ispec_is_elastic
  character(len=256) prname,LOCAL_PATH


!other parameters
  integer ::npoint_bound,ipoint,ipoint_temp,max_ndep_global,ndepth_table,&
            idepth,ndist_table,idist
  integer ::num_faces_this_mediatype
  integer ::size_depth_bound,size_dist_bound,max_ndist_iproc,max_ndist_allproc
  integer ::max_ndist_global
  integer, allocatable,dimension(:) ::ndist
  real(kind=CUSTOM_REAL) ::r_thisdepth,r_within_mesh
  integer ::izone,idiscontinuity,iglob_within_mesh
  integer,allocatable,dimension(:) ::depth_id,dist_id,iface_this_mediatype
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
  
  double precision ::theta_source,phi_source
  double precision,dimension(3) ::coord_source,coord_bound
  integer ::i,j,k,igll,iface,iface_temp,ispec,iglob

!parameters for sending and receiving
  integer::iproc,req_send,req_recv

!source
  integer :: yr,jda,ho,mi
  double precision :: lat_sour(1),long_sour(1),depth_sour(1)
  double precision ::t_cmt(1),hdur(1)
  double precision ::theta_sour,phi_sour
  double precision moment_tensor(6,1)

  double precision sec


  character(len=256) file_name,clean_LOCAL_PATH
  character(len=21) :: elastic_file
  character(len=22) :: acoustic_file
  data elastic_file /'id_depth_elastic.info'/
  data acoustic_file /'id_depth_acoustic.info'/

!  character(len=21) :: id_dep_ela_file
!  character(len=20) :: id_dist_ela_file
!  character(len=22) :: id_dep_acou_file
!  character(len=21) :: id_dist_acou_file
!  data id_dep_ela_file /'id_depth_elastic.info'/
!  data id_dist_ela_file /'id_dist_elastic.info'/
!  data id_dep_acou_file /'id_depth_acoustic.info'/
!  data id_dist_acou_file /'id_dist_elastic.info'/


!  print *,'read depth',media_type
!  open(unit=104,file=MF_IN_DATA_FILES_PATH(1:len_trim(MF_IN_DATA_FILES_PATH)) &
!       //'SEMtoTele_Par_file',status='old',action='read')
!  call read_value_double_precision(depth_tolerence,&
!       'SEMtoTele.depth_tolerence')
!  close(104)
!  if(myrank.eq.0) write(IMAIN,*) '  depth_torelence=',depth_tolerence


  num_faces_this_mediatype=0
  do iface=1,num_abs_boundary_faces
     ispec = abs_boundary_ispec(iface)
!     print *,'iface',iface,ispec_is_elastic(ispec)
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
  do iface=1,num_abs_boundary_faces
     ispec = abs_boundary_ispec(iface)
     if(ispec_is_acoustic(ispec).and.media_type.eq.2) then
        iface_temp=iface_temp+1
        iface_this_mediatype(iface_temp)=iface
        if(DEBUG_COUPLING)  print *,'acoustic ispec',iface
     else if(ispec_is_elastic(ispec).and.media_type.eq.1) then
        iface_temp=iface_temp+1
        iface_this_mediatype(iface_temp)=iface
     else if(.not.ispec_is_acoustic(ispec).and..not.ispec_is_elastic(ispec)) then
        print *,'media type',ispec_is_acoustic(ispec),ispec_is_elastic(ispec),media_type,&
               .not.ispec_is_acoustic(ispec).and..not.ispec_is_acoustic(ispec)
        call exit_MPI(myrank,'unknown or wrong media type')
     end if
  end do



  npoint_bound=NGLLX*NGLLY*num_faces_this_mediatype
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
  
!  if(DEBUG_COUPLING)  print *,'allocate finished',myrank


!to be fixed
!For COUP_TYPE=2 or 3, you need to uncomment the following line and make it
!correct
!  call get_cmt(yr,jda,ho,mi,sec,t_cmt,hdur,lat_sour,long_sour,depth_sour,moment_tensor,1)


  theta_sour = PI/2.0d0 - lat_sour(1)*PI/180.0d0
  phi_sour = long_sour(1)*PI/180.0d0
  coord_source(1)=dsin(theta_sour)*dcos(phi_sour)
  coord_source(2)=dsin(theta_sour)*dsin(phi_sour)
  coord_source(3)=dcos(theta_sour)

  ipoint=0
  do iface=1,num_faces_this_mediatype
     ispec = abs_boundary_ispec(iface_this_mediatype(iface))
     do igll = 1,NGLLSQUARE
       i = abs_boundary_ijk(1,igll,iface_this_mediatype(iface))
       j = abs_boundary_ijk(2,igll,iface_this_mediatype(iface))
       k = abs_boundary_ijk(3,igll,iface_this_mediatype(iface))
       iglob=ibool(i,j,k,ispec)
       x=xstore_dummy(iglob)
       y=ystore_dummy(iglob)
       z=zstore_dummy(iglob)
       r=dsqrt(x**2+y**2+z**2)
       ipoint=ipoint+1
       depth_bound(ipoint)=R_EARTH_SURF-r        
       coord_bound(1)=x/r;coord_bound(2)=y/r;coord_bound(3)=z/r

       dist_rad=dacos(dot_product(coord_bound,coord_source)*(1.d0-TINYVAL))
!WENBO
!       dist_bound(ipoint) = r*dist_rad
       dist_bound(ipoint) = dist_rad*180.d0/PI

       if(depth_bound(ipoint)<100.and.DEBUG_COUPLING) &
         print *,'check_dist_reslult',iface,igll,x,y,z,dist_rad*180.0/3.1415926,myrank
       if(depth_bound(ipoint)<5104900.0.and.myrank.eq.0.and.DEBUG_COUPLING) &
         print*,'check_topabs',iface,igll,i,j,k,dsqrt(x**2+y**2+z**2),dist_rad*180.0/3.1415926,myrank

     end do
  end do

max_ndep_global=MAX_NDEPTH_TABLE
!if(myrank.ne.0)  max_ndep_global=1
allocate(depth_table_global(max_ndep_global))

!construct depth table and get depth id
call find_global_id(depth_bound,size_depth_bound,npoint_bound,depth_tolerence,depth_id,&
                      depth_table_global,ndepth_table,max_ndep_global,&
                      NPROC,myrank)

ndiscon_found_deptable=0
if(ndepth_table.gt.0) then


 if(media_type.eq.1) then
     file_name = trim(LOCAL_PATH(1:len_trim(LOCAL_PATH)))//'/dist_table_elastic'
 else
     file_name = trim(LOCAL_PATH(1:len_trim(LOCAL_PATH)))//'/dist_table_acoustic'
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
        print *,myrank,"Error in counting coupoing points",ipoint_temp,npoint_bound
        call exit_MPI(myrank,"Error in counting coupoing points")
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
  if(izone_dep_table(1)-izone_dep_table(ndepth_table).ne.ndiscon_found_deptable) then
      print *,'check_izone',izone_dep_table(1),izone_dep_table(ndepth_table),&
            ndepth_table,ndiscon_found_deptable
      stop 'Error of discontinuity setting or depth table'
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
    if(idist.ne.ndist(idepth)) print *,myrank,"Error in counting DSM-SEM &
         coupling depths",idist,ndist(idepth)


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
  ipoint=0
  do iface=1,num_faces_this_mediatype
    do igll=1,NGLLSQUARE
     ipoint=ipoint+1
     if(idepth_is_discont(depth_id(ipoint)).eq.1) then
       ispec= abs_boundary_ispec(iface_this_mediatype(iface))
       i = abs_boundary_ijk(1,igll,iface_this_mediatype(iface))
       j = abs_boundary_ijk(2,igll,iface_this_mediatype(iface))
       k = abs_boundary_ijk(3,igll,iface_this_mediatype(iface))
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
         depth_id(ipoint)=depth_id(ipoint)+1
       end if
     else
       depth_id(ipoint)=depth_id(ipoint)+move_idepth(depth_id(ipoint))
     end if
    end do
  end do


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
    do igll=1,NGLLSQUARE
     ipoint=ipoint+1
     write(117,*) iface_this_mediatype(iface),igll,depth_id(ipoint),dist_id(ipoint),dist_bound(ipoint)
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
       file_name = trim(LOCAL_PATH(1:len_trim(LOCAL_PATH)))//'/depth_table_elastic'
     else
       file_name = trim(LOCAL_PATH(1:len_trim(LOCAL_PATH)))//'/depth_table_acoustic'
     end if
     open(17,file=file_name(1:len_trim(file_name)),status='unknown',action='write',form='formatted')
     write(17,*) NGLLZ
     write(17,*) wzgll(1:NGLLZ)/2.d0
     write(17,*) ndepth_table+ndiscon_found_deptable
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



subroutine find_global_id(input_array,size_array,npoints_array,tolerence_value,id_inglobal,&
                       table_global,n_table_global,max_ntable_global,NPROC,myrank)
use constants
implicit none
!include "constants.h"

!input and output
integer,intent(in) ::NPROC,myrank
integer,intent(in) ::size_array,npoints_array,max_ntable_global
real(kind=CUSTOM_REAL),dimension(size_array),intent(in) ::input_array
double precision, intent(in) ::tolerence_value

integer,intent(out) ::n_table_global
integer,dimension(size_array),intent(out) ::id_inglobal
real(kind=CUSTOM_REAL),dimension(max_ntable_global),intent(out) ::table_global

!local parameters
integer ::n_table,max_n_table,ntable_all_local,ipoint,itable_test,itable,itable_temp
integer ::index_iproc,iproc
integer, dimension(:),allocatable ::id_inlocal
real(kind=CUSTOM_REAL),dimension(:),allocatable ::table_local,table_all_local
integer, dimension(:), allocatable ::map_localtoglobal,map_all
integer,dimension(:), allocatable ::ntable_eachproc
integer ::req_send,req_recv
integer,dimension(1) ::ntable_buffer,n_table_global_buffer

n_table_global=0
id_inglobal(:)=0
table_global(:)=0

max_n_table=npoints_array
if(npoints_array.gt.0) then
   allocate(table_local(max_n_table))
   allocate(id_inlocal(npoints_array))
end if
if(myrank.eq.0) then
   allocate(ntable_eachproc(NPROC))
   ntable_eachproc(:)=0
end if
  n_table=0
  do ipoint=1,npoints_array
      itable_test=1
      do while(itable_test.le.n_table.and.&
               abs(table_local(itable_test)-input_array(ipoint))&
               .gt.tolerence_value)
                   itable_test=itable_test+1
      end do
!     if(DEBUG_COUPLING) print *,'test',myrank,ipoint,idepth_test,ndepth_table,&
!             depth_table_local(idepth_test)-depth_bound(ipoint),depth_tolerence
      if(itable_test.gt.n_table) then
           n_table=itable_test
           table_local(n_table)=input_array(ipoint)
!           print *,'local_table',depth_table_local(ndepth_table),ndepth_table
      end if
      id_inlocal(ipoint)=itable_test
!     if(DEBUG_COUPLING) print *,'id_inlocal',id_inlocal(ipoint),depth_bound(ipoint),myrank,ipoint
  end do
  if(DEBUG_COUPLING) print *,'ndepth_table1',n_table,table_local(:),myrank
  if(n_table.gt.0)  allocate(map_localtoglobal(n_table))

!send and receive ndepth_table
  if(myrank.eq.0) then
!        if(DEBUG_COUPLING) print *,'recv0'
         ntable_eachproc(1)=n_table
         ntable_all_local=ntable_eachproc(1)
    do iproc=1,NPROC-1
!         if(DEBUG_COUPLING) print *,'recv ndepth',iproc
!         call recv_i(ntable_eachproc(iproc+1),1,iproc,0)
         call recv_i(ntable_buffer,1,iproc,0)
         ntable_eachproc(iproc+1)=ntable_buffer(1)
!        if(DEBUG_COUPLING) print *,'recv ndepth1',iproc,ndepth_eachproc(iproc+1)
         ntable_all_local=ntable_all_local+ntable_eachproc(iproc+1)
    end do
    if(ntable_all_local.gt.0) allocate(table_all_local(ntable_all_local))
  else
!    call send_i((/n_table/),1,0,0)
    ntable_buffer(1)=n_table
!    call send_i(n_table,1,0,0)
    call send_i(ntable_buffer,1,0,0)
  end if
!send and receive depth_table_temp
  if(myrank.eq.0) then
!copy local table of rank0
        if(n_table.ge.1) table_all_local(1:n_table)=table_local(1:n_table)
         index_iproc=n_table+1
!receiver other ranks' local table
    do iproc=1,NPROC-1
         if(ntable_eachproc(iproc+1).gt.0) then
           req_recv=iproc
           call irecv_cr(table_all_local(index_iproc),ntable_eachproc(iproc+1),&
                         iproc,1,req_recv)
           call wait_req(req_recv)
           index_iproc=index_iproc+ntable_eachproc(iproc+1)
         end if
    end do
  else
    req_send=myrank
    if(n_table.gt.0) then
!           if(DEBUG_COUPLING)  print *,'n_table',n_table,max_n_table
            call isend_cr(table_local,n_table,&
                                        0,1,req_send)
            call wait_req(req_send)
    end if
  end if
!merge the local depth tables, and then construct map_localtoglobal
  if(myrank.eq.0) then
!   if(DEBUG_COUPLING) print *,'depth_all',table_all_local(1:ntable_all_local)
    do itable=1,ntable_all_local
      if(abs(table_all_local(itable)).lt.tolerence_value) then
!       if(DEBUG_COUPLING) print *,'depth0',table_all_local(itable)
!WENBO
!        table_all_local(itable)=tolerence_value
        table_all_local(itable)=0.0
      else if(table_all_local(itable).lt.-tolerence_value) then
        call exit_MPI(myrank,'Error of input array!!! array(i)<-tolerence_value. &
                              Usually, it is negative depths.')
      end if
    end do
!    max_ntable_global=ntable_all_local
    if(max_ntable_global.gt.0) then
!      allocate(table_global(max_ntable_global))
      allocate(map_all(ntable_all_local))
    end if
    table_global(:)=1.e25
    map_all(:)=0
    do itable=1,ntable_all_local
      itable_test=1
      do while(itable_test.le.n_table_global.and.&
               table_global(itable_test).lt.table_all_local(itable) )
                   itable_test=itable_test+1
      end do
      if(itable_test.eq.1) then
        if(abs(table_global(itable_test)-table_all_local(itable)).lt. &
               tolerence_value) then
         map_all(itable)=itable_test
        else
         do itable_temp=n_table_global,itable_test,-1
          if(itable_temp.lt.-1) stop 'Error itable_temp'
          table_global(itable_temp+1)=table_global(itable_temp)
         end do
         do itable_temp=1,itable-1
          if(map_all(itable_temp).ge.itable_test) then
             map_all(itable_temp)=map_all(itable_temp)+1
          end if
         end do
         table_global(itable_test)=table_all_local(itable)
         map_all(itable)=itable_test
         n_table_global=n_table_global+1
         if(n_table_global.gt.max_ntable_global) &
             call exit_MPI(myrank,'Error,n_table_global>max_ntable_global')
          if(DEBUG_COUPLING) print *,'case3',table_global(itable_test),table_all_local(itable),&
             itable_test,n_table_global
!        if(DEBUG_COUPLING) print *,'map',map_all(:)

        end if
      else
        if(abs(table_global(itable_test)-table_all_local(itable)).lt. &
               tolerence_value) then
         map_all(itable)=itable_test
         if(DEBUG_COUPLING) print *,'case1',table_global(itable_test),table_all_local(itable)
        else if(itable_test.gt.1.and.itable_test.le.n_table_global+1.and. &
             (abs(table_global(itable_test-1)-table_all_local(itable)) &
                .lt.tolerence_value) ) then
          map_all(itable)=itable_test-1
         if(DEBUG_COUPLING) print *,'case2',table_global(itable_test),table_all_local(itable)
        else
         do itable_temp=n_table_global,itable_test,-1
          if(itable_temp.lt.-1) stop 'Error itable_temp' 
          table_global(itable_temp+1)=table_global(itable_temp)
         end do
         do itable_temp=1,itable-1
          if(map_all(itable_temp).ge.itable_test) then
             map_all(itable_temp)=map_all(itable_temp)+1
          end if
         end do
         table_global(itable_test)=table_all_local(itable)
         map_all(itable)=itable_test
         n_table_global=n_table_global+1
         if(n_table_global.gt.max_ntable_global) &
             call exit_MPI(myrank,'Error,n_table_global>max_ntable_global')
          if(DEBUG_COUPLING) print *,'case3',table_global(itable_test),table_all_local(itable),&
             itable_test,n_table_global
!        if(DEBUG_COUPLING) print *,'map',map_all(:)

        end if
       end if
!      if(DEBUG_COUPLING) print *,'map_all',idep,map_all(idep),idepth_test,ndepth_table_global,&
!               depth_all_local(idep)
    end do
  end if
!call bcast_all_i((/n_table_global/),1)
!call bcast_all_i(n_table_global,1)
n_table_global_buffer(1)=n_table_global
call bcast_all_i(n_table_global_buffer,1)
n_table_global=n_table_global_buffer(1)
!send map_localtoglobal
  if(myrank.eq.0) then
!    if(DEBUG_COUPLING) print *,'table',map_all(:)
         index_iproc=n_table+1
    do iproc=1,NPROC-1
        if(ntable_eachproc(iproc+1).gt.0) then
          req_send=iproc
          call isend_i(map_all(index_iproc),ntable_eachproc(iproc+1), &
               iproc,2,req_send)
          call wait_req(req_send)
          index_iproc=index_iproc+ntable_eachproc(iproc+1)
        end if
    end do
    if(n_table.gt.0) map_localtoglobal(1:n_table)=map_all(1:n_table)
  else
    if(n_table.gt.0) then
!        if(myrank.eq.61.and.DEBUG_COUPLING) print *,'ndepth_table',ndepth_table,map_localtoglobal
        req_recv=myrank
        call irecv_i(map_localtoglobal,n_table,0,2,req_recv)
        call wait_req(req_recv)
    end if
  end if

!construct global id table
  do ipoint=1,npoints_array
    id_inglobal(ipoint)= map_localtoglobal(id_inlocal(ipoint))
!    if(id_inglobal(ipoint).eq.0.and.DEBUG_COUPLING) &
!      print *,'id_depth',ipoint,myrank,id_inlocal(ipoint),n_table,npoints_array,map_localtoglobal(:)
  end do

if(npoints_array.gt.0) then
     deallocate(table_local)
     deallocate(id_inlocal)
end if

if(myrank.eq.0) then
   deallocate(ntable_eachproc)
   if(ntable_all_local.gt.0) deallocate(table_all_local)
   if(max_ntable_global.gt.0)  deallocate(map_all)
end if

if(n_table.gt.0)  deallocate(map_localtoglobal)


end subroutine
