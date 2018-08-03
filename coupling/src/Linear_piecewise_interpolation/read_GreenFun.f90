subroutine read_GreenFunc_elastic(ifreq)
    use convolution_par
    use coupling_SEM_DSM_par, only:myrank,do_coupling_RTZ
    use constants
    implicit none

! other variables
    character(len=80) ::green_disp_file,green_epsilon_file
    integer ::icomp,icomp_force,ifreq,ios
    character(len=80),dimension(ncomp) ::force_comp
    data force_comp /'zcomp','rcomp','tcomp'/


   do icomp_force=1,ncomp
    if(do_coupling_RTZ(icomp_force)) then
      write(green_disp_file,"('./DSM_input/',a5,'/disp_solid/freq_',i5.5)") force_comp(icomp_force),ifreq-1
!     open(unit=15,file=trim(green_disp_file),status='unknown', form='formatted',iostat=ios)
      open(unit=15,file=trim(green_disp_file),status='unknown', form='unformatted',iostat=ios)
      if(ios /= 0) stop 'error reading disp_solid'
!code- single force
!     read(15)disp_Green(:,:,icomp_force)

!code-velo_stress_quadra
      do icomp=1,3
        read(15) disp_Green_elas(icomp,:,icomp_force)
      end do
      close(15)
    end if
   end do

   do icomp_force=1,3
    if(do_coupling_RTZ(icomp_force)) then
       write(green_epsilon_file,"('./DSM_input/',a5,'/sigma/freq_',i5.5)") force_comp(icomp_force),ifreq-1
!      open(unit=16,file=trim(green_epsilon_file),status='unknown',form='formatted',iostat=ios)
       open(unit=16,file=trim(green_epsilon_file),status='unknown',form='unformatted',iostat=ios)
       if(ios /= 0) stop 'error reading GREENS_EPSILON'
       do icomp=1,ncomp_stress
          read(16)epsilon_Green(:,icomp,icomp_force)
       end do
       close(16)
     end if
    end do


end subroutine read_GreenFunc_elastic

subroutine read_GreenFunc_acoustic(ifreq)
    use convolution_par
    use coupling_SEM_DSM_par, only:myrank,do_coupling_RTZ
    use constants
    implicit none

! other variables
    character(len=80) ::green_disp_file,green_pressure_file
    integer ::icomp,icomp_force,ifreq,ios
    character(len=80),dimension(ncomp) ::force_comp
    data force_comp /'zcomp','rcomp','tcomp'/


   do icomp_force=1,ncomp
    if(do_coupling_RTZ(icomp_force)) then
      write(green_disp_file,"('./DSM_input/',a5,'/disp_fluid/freq_',i5.5)") force_comp(icomp_force),ifreq-1
!     open(unit=15,file=trim(green_disp_file),status='unknown', form='formatted',iostat=ios)
      open(unit=15,file=trim(green_disp_file),status='unknown', form='unformatted',iostat=ios)
      if(ios /= 0) stop 'error reading disp_fluid'
!code- single force
!     read(15)disp_Green(:,:,icomp_force)

!code-velo_stress_quadra
      do icomp=1,3
        read(15) disp_Green_acous(icomp,:,icomp_force)
      end do
      close(15)
    end if
   end do

   do icomp_force=1,ncomp
    if(do_coupling_RTZ(icomp_force)) then
      write(green_pressure_file,"('./DSM_input/',a5,'/pressure/freq_',i5.5)") force_comp(icomp_force),ifreq-1
!     open(unit=16,file=trim(green_pressure_file),status='unknown',form='formatted',iostat=ios)
      open(unit=16,file=trim(green_pressure_file),status='unknown',form='unformatted',iostat=ios)
      if(ios /= 0) stop 'error reading GREENS_EPSILON'
      read(16)pressure_Green(:,icomp_force)
      close(16)
    end if
  end do

end subroutine read_GreenFunc_acoustic
