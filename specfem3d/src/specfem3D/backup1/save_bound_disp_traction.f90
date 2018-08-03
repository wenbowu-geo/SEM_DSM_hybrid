subroutine save_bound_disp_traction()
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use SEMtoTele_par
  implicit none

  integer ::iface_bound,iface_bound_elastic
  integer ::ipoint_elastic,ipoint_elas_acous,ipoint_iface
  integer ::imax,imin,di,jmax,jmin,dj,kmax,kmin,dk
  integer ::i,j,k,it_bound,it_bound_tmp,ispec,nstep_iblock
  integer ::ipackage
  integer ::ier1,ier2
  character(len=256) ::final_LOCAL_PATH,clean_LOCAL_PATH
  character(len=256) ::dir_this_proc,disp_bound_file,&
                       traction_bound_file
  
  it_bound=(it-1)/DECIMATE_COUPLING+1

!  if( USE_OUTPUT_FILES_PATH ) then
!      final_LOCAL_PATH = OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH))! //'/'
!  else
      ! suppress white spaces if any
      clean_LOCAL_PATH = adjustl(LOCAL_PATH)
      ! create full final local path
      final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH))! //'/'
!  endif
  dir_this_proc=trim(final_LOCAL_PATH)
  write(dir_this_proc,"(a,'iproc',i4.4)")trim(dir_this_proc),myrank

  ipackage=0
  iface_bound_elastic=0
  do iface_bound=1,Nele_Bound

   ispec=Bound_Info(iface_Bound)%ispec_bound
   imin=Bound_Info(iface_Bound)%istart;imax=Bound_Info(iface_Bound)%iend
   jmin=Bound_Info(iface_Bound)%jstart;jmax=Bound_Info(iface_Bound)%jend
   kmin=Bound_Info(iface_Bound)%kstart;kmax=Bound_Info(iface_Bound)%kend

   if(ispec_is_elastic(ispec)) iface_bound_elastic=iface_bound_elastic+1
   di=1;dj=1;dk=1

   ipoint_iface=0
   do k=kmin,kmax,dk
     do j=jmin,jmax,dj
       do i=imin,imax,di
!         if(ispec_is_elastic(ispec)) then
            if(LOW_RESOLUTION) then
                   ipoint_elastic=iface_bound_elastic
                   ipoint_elas_acous=iface_bound
            else
                   ipoint_iface=ipoint_iface+1
                   ipoint_elas_acous=(iface_bound-1)*NGLLX*NGLLY+ipoint_iface
                   ipoint_elastic=(iface_bound_elastic-1)*NGLLX*NGLLY+ipoint_iface
            end if

            if(mod(ipoint_elas_acous,NPOINTS_PER_PACK).eq.1) then
               ipackage=ipackage+1
               write(disp_bound_file,"(a,'/disp_pack',i6.6)")trim(dir_this_proc),ipackage
               write(traction_bound_file,"(a,'/traction_pack',i6.6)")& 
                                                             trim(dir_this_proc),ipackage

               if(it_bound.le.NSTEP_BETWEEN_OUTPUTBOUND) then
!                 open(unit=11,file=trim(disp_bound_file),status='unknown',&
!                   action='write',form='formatted',iostat=ier1)
!                 open(unit=21,file=trim(traction_bound_file),status='unknown',&
!                   action='write',form='formatted',iostat=ier2)
                 open(unit=11,file=trim(disp_bound_file),status='unknown',&
                   action='write',access='stream',form='unformatted',iostat=ier1)
                 open(unit=21,file=trim(traction_bound_file),status='unknown',&
                   action='write',access='stream',form='unformatted',iostat=ier2)
               else 
!                 open(unit=11,file=trim(disp_bound_file),status='old',&
!                   action='write',form='formatted',position='append',iostat=ier1)
!                 open(unit=21,file=trim(traction_bound_file),status='old',&
!                   action='write',form='formatted',position='append',iostat=ier2)
                 open(unit=11,file=trim(disp_bound_file),status='old',&
                   action='write',access='stream',form='unformatted',position='append',iostat=ier1)
                 open(unit=21,file=trim(traction_bound_file),status='old',&
                   action='write',access='stream',form='unformatted',position='append',iostat=ier2)
               end if
       
               if( ier1 /= 0 ) then
                 print*,'error: could not open  file'
                 print*,'path:',disp_bound_file(1:len_trim(disp_bound_file))
                 call exit_mpi(myrank,'error opening disp_bound file')
               endif

               if( ier2 /= 0 ) then
                 print*,'error: could not open  file'
                 print*,'path:',traction_bound_file(1:len_trim(traction_bound_file))
                 call exit_mpi(myrank,'error opening traction_bound file')
               endif
             end if
            if(it.eq.NSTEP) then
              if(DECIMATE_COUPLING.eq.1) then
                nstep_iblock=mod(it_bound,NSTEP_BETWEEN_OUTPUTBOUND)
                if(nstep_iblock.eq.0) nstep_iblock=NSTEP_BETWEEN_OUTPUTBOUND
              else 
                if(mod(it_bound,NSTEP_BETWEEN_OUTPUTBOUND).ne.0) then
                  nstep_iblock=mod(it_bound,NSTEP_BETWEEN_OUTPUTBOUND)
                else if(it-1.eq.(it_bound-1)*DECIMATE_COUPLING) then
                  nstep_iblock=NSTEP_BETWEEN_OUTPUTBOUND
                else if(it-((it_bound-1)*DECIMATE_COUPLING+1).lt.DECIMATE_COUPLING) then
                  nstep_iblock=0
                else
                    call exit_MPI(myrank,'Error of nstep_iblock!')
                end if
              end if
            else
                nstep_iblock=NSTEP_BETWEEN_OUTPUTBOUND
            end if
            if(myrank.eq.0.and.ipoint_elas_acous.eq.1) &
                    print *,'disp1_bound',it,it_bound,nstep_iblock,ipackage
            do it_bound_tmp=1,nstep_iblock
!               write(11,*) disp_bound(it_bound_tmp,:,ipoint_elas_acous)
!               write(21,*) traction_bound(it_bound_tmp,:,ipoint_elas_acous)
               write(11) disp_bound(it_bound_tmp,:,ipoint_elas_acous)
               write(21) traction_bound(it_bound_tmp,:,ipoint_elas_acous)
            end do
            if(mod(ipoint_elas_acous,NPOINTS_PER_PACK).eq.0.or.ipoint_elas_acous.eq.npoints_bound) then
               close(11)
               close(21)
            end if
!        end if
       end do
     end do
   end do

  end do


end subroutine save_bound_disp_traction
