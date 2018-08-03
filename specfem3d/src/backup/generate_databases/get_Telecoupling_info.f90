subroutine get_Telecoupling_info(myrank,xstore_dummy,ystore_dummy,zstore_dummy,ibool,nspec,&
            nglob,kappastore,mustore,LOCAL_PATH,&
            c11store,c12store,c13store,c14store,c15store,c16store,&
            c22store,c23store,c24store,c25store,c26store,c33store,&
            c34store,c35store,c36store,c44store,c45store,c46store,&
            c55store,c56store,c66store,ANISOTROPY,NSPEC_ANISO,&
            xigll,yigll,zigll,wxgll,wygll,wzgll,&
            wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,ispec_is_elastic,&
            ispec_is_acoustic)
  use SEMtoTele_par
  use generate_databases_par, only:NGNOD2D,NGNOD

  implicit none

!input parameters
integer ::nspec,nglob,NSPEC_ANISO,myrank
character(len=256) LOCAL_PATH
logical ::ANISOTROPY
logical, dimension(nspec) :: ispec_is_elastic,ispec_is_acoustic
double precision, dimension(NGLLX,NGLLY) :: wgllwgll_xy
double precision, dimension(NGLLX,NGLLZ) :: wgllwgll_xz
double precision, dimension(NGLLY,NGLLZ) :: wgllwgll_yz
real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec)::kappastore,mustore
real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO) :: &
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store, &
            c34store,c35store,c36store,c44store,c45store,c46store, &
            c55store,c56store,c66store
integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy
double precision, dimension(NGLLX)::xigll,wxgll
double precision, dimension(NGLLY)::yigll,wygll
double precision, dimension(NGLLZ)::zigll,wzgll




!other variables
integer ::ispec,iface,ipoint_elas_acous,ipoint_elastic
integer ::ipoint_iface,iface_bound,iface_bound_elastic
integer,target ::i,j,k
integer ::ipackage
integer        ::face_type
integer ::imax,imin,di,jmax,jmin,dj,kmax,kmin,dk
integer, pointer::index1,index2


double precision, dimension(NDIM2D,NGNOD2D,NGLLY,NGLLZ):: dershape2D_x
double precision, dimension(NDIM2D,NGNOD2D,NGLLX,NGLLZ):: dershape2D_y
double precision, dimension(NDIM2D,NGNOD2D,NGLLX,NGLLY):: dershape2D_bottom
double precision, dimension(NDIM2D,NGNOD2D,NGLLX,NGLLY):: dershape2D_top
double precision, dimension(NGNOD2D,NGLLY,NGLLZ)       :: shape2D_x,shape2D_y,&
                                                          shape2D_bottom,shape2D_top
!double precision, dimension(NGNOD2D)                   ::xcoord,ycoord,zcoord
real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY)          :: jacobian2Dw_face
real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ)    :: jacobian2D
real(kind=CUSTOM_REAL)                                 :: jacobian_whole_face
real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY)     ::normal_face

real(kind=CUSTOM_REAL) ::lambdal,mu
real(kind=CUSTOM_REAL) ::c11,c12,c13,c14,c15,c16,&
            c22,c23,c24,c25,c26,c33,&
            c34,c35,c36,c44,c45,c46,&
            c55,c56,c66


 character(len=256) ::normal_file_solid,jacobian_file_solid,elas_coef_file_solid,&
                      coord_file_solid
 character(len=256) ::boundinfo_file_solid
 character(len=256) ::final_LOCAL_PATH,clean_LOCAL_PATH
 character(len=256) ::dir_this_proc
 integer ::ier


! if( USE_OUTPUT_FILES_PATH ) then
!      final_LOCAL_PATH = OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH))
! else
      ! suppress white spaces if any
      clean_LOCAL_PATH = adjustl(LOCAL_PATH)
      ! create full final local path
      final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) 
! endif
 dir_this_proc=trim(final_LOCAL_PATH)
 write(dir_this_proc,"(a,'iproc',i4.4,'/')")trim(dir_this_proc),myrank

 call system('mkdir -p '// adjustl(trim( dir_this_proc ) ) )


! get the 2-D shape functions
call get_shape2D(myrank,shape2D_x,dershape2D_x,yigll,zigll,NGLLY,NGLLZ,NGNOD,NGNOD2D)
call get_shape2D(myrank,shape2D_y,dershape2D_y,xigll,zigll,NGLLX,NGLLZ,NGNOD,NGNOD2D)
call get_shape2D(myrank,shape2D_bottom,dershape2D_bottom,xigll,yigll,NGLLX,NGLLY,NGNOD,NGNOD2D)
call get_shape2D(myrank,shape2D_top,dershape2D_top,xigll,yigll,NGLLX,NGLLY,NGNOD,NGNOD2D)


ipackage=0
iface_bound_elastic=0

do iface_bound=1,Nele_Bound
    face_type=Bound_Info(iface_bound)%face_type
    select case(face_type)
    case(1)
      index1=>j;index2=>k
      iface=1
    case(2)
      index1=>j;index2=>k
      iface=2
    case(3)
      index1=>i;index2=>k
      iface=3
    case(4)
      index1=>i;index2=>k
      iface=4
   case(5)
     index1=>i;index2=>j
     iface=5
   case(6)
     index1=>i;index2=>j
     iface=6
   case default
     stop 'Error boundary type'
   end select


   ispec=Bound_Info(iface_bound)%ispec_bound
   if(ispec_is_elastic(ispec)) iface_bound_elastic=iface_bound_elastic+1

   call get_jacobian_boundary_face(myrank,nspec, &
               xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
               dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top,&
               dble(wgllwgll_xy),dble(wgllwgll_xz),dble(wgllwgll_yz),&
               ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)
!   do i=1,NGLLX
!      do j=1,NGLLY
!         call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
!                                      ibool,nspec,nglob, &
!                                      xstore,ystore,zstore, &
!                                      normal_face(:,i,j) )
!      end do
!   end do


   if(LOW_RESOLUTION) then
     jacobian_whole_face=0.d0
     do i=1,NGLLX
       do j=1,NGLLY
           jacobian_whole_face=jacobian_whole_face+jacobian2Dw_face(i,j)
           if(myrank.eq.0.and.iface_bound.eq.1) print *,'jacob_whole',i,j,k,&
                 jacobian2Dw_face(i,j),jacobian_whole_face 
       end do
     end do
   end if
   imin=Bound_Info(iface_bound)%istart;imax=Bound_Info(iface_bound)%iend
   jmin=Bound_Info(iface_bound)%jstart;jmax=Bound_Info(iface_bound)%jend
   kmin=Bound_Info(iface_bound)%kstart;kmax=Bound_Info(iface_bound)%kend
   di=1;dj=1;dk=1
   
   ipoint_iface=0
   do k=kmin,kmax,dk
     do j=jmin,jmax,dj
       do i=imin,imax,di
         if(ispec_is_elastic(ispec)) then
            if(.not.ANISOTROPY) then
              mu=mustore(i,j,k,ispec)
              lambdal=kappastore(i,j,k,ispec)-2.d0/3.d0*mustore(i,j,k,ispec)
              c11=lambdal+TWO*mu; c12=lambdal; c13=lambdal
              c14=0.d0;           c15=0.d0;   c16=0.d0
              c22=lambdal+TWO*mu; c23=lambdal; c24=0.d0
              c25=0.d0;           c26=0.d0;   
              c33=lambdal+TWO*mu; c34=0.d0;    c35=0.d0
              c36=0.d0
              c44=mu;             c45=0.d0;    c46=0.d0
              c55=mu;             c56=0.d0
              c66=mu
            else
              c11=c11store(i,j,k,ispec); c12=c12store(i,j,k,ispec)
              c13=c13store(i,j,k,ispec); c14=c14store(i,j,k,ispec)
              c15=c15store(i,j,k,ispec); c16=c16store(i,j,k,ispec)
              c22=c22store(i,j,k,ispec); c23=c23store(i,j,k,ispec)
              c24=c24store(i,j,k,ispec); c25=c25store(i,j,k,ispec)
              c26=c26store(i,j,k,ispec); c33=c33store(i,j,k,ispec)
              c34=c34store(i,j,k,ispec); c35=c35store(i,j,k,ispec)
              c36=c36store(i,j,k,ispec); c44=c44store(i,j,k,ispec)
              c45=c45store(i,j,k,ispec); c46=c46store(i,j,k,ispec)
              c55=c55store(i,j,k,ispec); c56=c56store(i,j,k,ispec)
              c66=c66store(i,j,k,ispec)
            end if
         end if

           if(LOW_RESOLUTION) then
               ipoint_elas_acous=iface_bound
               ipoint_elastic=iface_bound_elastic
           else
               ipoint_iface=ipoint_iface+1
               ipoint_elas_acous=(iface_bound-1)*NGLLX*NGLLY+ipoint_iface
               ipoint_elastic=(iface_bound_elastic-1)*NGLLX*NGLLY+ipoint_iface
           end if
!           print *,'normal',iface,normal_face(:,index1,index2),&
!                   xstore(ibool(i,j,k,ispec)),ystore(ibool(i,j,k,ispec)),&
!                   zstore(ibool(i,j,k,ispec))
           if(iface.eq.5.or.iface.eq.1.or.iface.eq.4) then
              normal_face(:,index1,index2)=-normal_face(:,index1,index2)
           end if
           normal_vect_bound(:,ipoint_elas_acous)=normal_face(:,index1,index2)

           if(mod(ipoint_elas_acous,NPOINTS_PER_PACK).eq.1) then
               ipackage=ipackage+1
               write(boundinfo_file_solid,"(a,'bound_info_pack',i6.6)")trim(dir_this_proc),ipackage
               open(unit=115,file=trim(boundinfo_file_solid),status='unknown',&
                   action='write',form='formatted',iostat=ier)
               if( ier /= 0 ) then
                 print*,'error: could not open  file'
                 print*,'path:',boundinfo_file_solid(1:len_trim(boundinfo_file_solid))
                 call exit_mpi(myrank,'error opening boundinfo_file_solid file')
               endif

           end if

!           write(115,*)id_depth_bound(ipoint),id_dist_bound(ipoint),ispec_is_elastic(ispec),ispec_is_acoustic(ispec)
           write(115,*)id_depth_bound(ipoint_elas_acous),ispec_is_elastic(ispec),ispec_is_acoustic(ispec)
!           write(115,*)c11,c12,c13,c14,c15,c16,&
!             c22,c23,c24,c25,c26,c33,&
!             c34,c35,c36,c44,c45,c46,&
!             c55,c56,c66
           write(115,*) normal_face(:,index1,index2)
           if(LOW_RESOLUTION) then
             write(115,*) jacobian_whole_face
           else
             write(115,*) jacobian2Dw_face(index1,index2)
           end if
           write(115,*) xstore_dummy(ibool(i,j,k,ispec)),ystore_dummy(ibool(i,j,k,ispec)),&
                       zstore_dummy(ibool(i,j,k,ispec))
!                       6371000.0-dsqrt(xstore_dummy(ibool(i,j,k,ispec))**2+&
!                       ystore_dummy(ibool(i,j,k,ispec))**2+zstore_dummy(ibool(i,j,k,ispec))**2)
           if(mod(ipoint_elas_acous,NPOINTS_PER_PACK).eq.0.or.ipoint_elas_acous.eq.npoints_bound) then
               close(115)
           end if

       end do
     end do
   end do
end do


end subroutine get_Telecoupling_info
