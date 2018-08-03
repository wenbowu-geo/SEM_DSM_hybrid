subroutine SEMtoTele() 
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use SEMtoTele_par
  implicit none
  integer ::it_bound

  if(it.eq.1) then
     call read_EleBound()
     print *,'read TelePara'
     call read_parameter_SEMtoTele(myrank)
     call compute_pml_rotation_matrix(myrank)
!     call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,&
!          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)
     print *,'set Varia'
     call set_Variables_Bound()
!     call save_package_list()

!get_Telecoupling_info is executed in xgene
     call get_Telecoupling_info()
  end if
!  print *,'calculate_strain_ela start',myrank,it
  if(DECIMATE_COUPLING.eq.1.or.mod(it,DECIMATE_COUPLING).eq.1) then
     it_bound=(it-1)/DECIMATE_COUPLING+1
     call compute_traction_disp_elastic()
  if(myrank.eq.40) print *,'calculate_strain_ela',myrank,it
     call compute_pres_disp_acoustic()
!  call calculate_strain_acoustic()
!  print *,'calculate_strain_acou',myrank,it
     !call compute_pres_disp_acoustic()
     if(it.eq.1.and.myrank.eq.195) print *,'coord_disp1',&
                      xstore(ibool(2,2,1,604)),ystore(ibool(2,2,1,604)),&
                      zstore(ibool(2,2,1,604))
     if(it.eq.1.and.myrank.eq.195) print *,'coord_disp2',&
                      xstore(ibool(3,3,1,604)),ystore(ibool(3,3,1,604)),&
                      zstore(ibool(3,3,1,604))

     if(myrank.eq.195) print *,'bound_disp1',displ(:,ibool(2,2,1,604))
     if(myrank.eq.195) print *,'bound_disp2',displ(:,ibool(3,3,1,604))
   else
     it_bound=-1
   end if
!  print *,'velocity',myrank,it
   if( (it.gt.1.and.mod(it_bound,NSTEP_BETWEEN_OUTPUTBOUND).eq.0 )&
      .or. it .eq. NSTEP) then
!WENBO 
!not comment it!
     if(myrank.eq.40) print *,"save_bound",it
     call  save_bound_disp_traction()
   end if

end subroutine SEMtoTele
