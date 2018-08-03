!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! generic tomography file
!
! note: the idea is to use an external, tomography velocity model
!
! most of the routines here are place-holders, please add/implement your own routines
!

  module tomography
  use constants
  !include "constants.h"

  ! for external tomography:
  ! file must be in ../in_data/files/ directory
  ! (regular spaced, xyz-block file in ascii)
  !character (len=80) :: TOMO_FILENAME = 'veryfast_tomography_abruzzo_complete.xyz'
  character (len=80) :: TOMO_FILENAME = 'tomography_model.xyz'
  character (len=80) :: DSM_FILENAME = 'dsm_input_file'

  !DSM input model
  integer ::n_structure_zone
  integer, dimension(:),allocatable::fluid_thiszone
  real(kind=CUSTOM_REAL),dimension(:),allocatable ::rmin_structure_zone,rmax_structure_zone,qmu_structure_zone,qkappa_structure_zone
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable ::rho_structure_zone,&
               vpv_structure_zone,vph_structure_zone,vsv_structure_zone,&
               vsh_structure_zone,eta_structure_zone




  ! model dimensions
  double precision :: ORIG_X,ORIG_Y,ORIG_Z
  double precision :: END_X,END_Y,END_Z
  double precision :: SPACING_X,SPACING_Y,SPACING_Z

  ! model parameter records
  real(kind=CUSTOM_REAL), dimension (:), allocatable :: vp_tomography,vs_tomography,rho_tomography,z_tomography

  !topography
  double precision ::ORIG_LON,ORIG_LAT, END_LON, END_LAT
  double precision ::SPACING_LON, SPACING_LAT
  integer :: NLON, NLAT
  real(kind=CUSTOM_REAL), dimension (:,:), allocatable :: topography

  !trench
  integer ::NTRENCH
  double precision,dimension(:),allocatable :: lat_trench,lon_trench

  ! model entries
  integer :: NX,NY,NZ
  integer :: nrecord

  ! min/max statistics
  double precision :: VP_MIN,VS_MIN,RHO_MIN,VP_MAX,VS_MAX,RHO_MAX

  end module tomography

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_tomography_broadcast(myrank)

  implicit none

  ! include "constants.h"
  ! include "precision.h"
  ! include 'mpif.h'
  integer :: myrank

  ! all processes read in same file
  ! note: for a high number of processes this might lead to a bottleneck
  call read_model_tomography(myrank)

  ! otherwise:

  ! only master reads in model file
  !if(myrank == 0) call read_external_model()
  ! broadcast the information read on the master to the nodes, e.g.
  !call MPI_BCAST(nrecord,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !if( myrank /= 0 ) allocate( vp_tomography(1:nrecord) )
  !call MPI_BCAST(vp_tomography,size(vp_tomography),CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  end subroutine model_tomography_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_tomography(myrank)

! start magnoni 29/11/09
! read Vp Vs and rho from extracted text file

! assuming that only tomography undefined material is allowed....
! and all the tomographic regions are collect inside one file called TOMO_FILENAME with homogenous resolution
! this could be problematic for example if the tomographic regions have different resolution
! leading to a waste of memory and cpu time in the partitioning process

  use tomography

  implicit none

  integer :: myrank

  ! local parameters
!  real(kind=CUSTOM_REAL) :: x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: array_read_tmp
  integer :: irecord,ier
  integer ::iz,iy,irecord0
  real(kind=CUSTOM_REAL) ::z_tmp
  character(len=256):: filename

  !TOMO_FILENAME='DATA/veryfast_tomography_abruzzo_complete.xyz'
  ! probably the simple position for the filename is the constat.h
  ! but it is also possible to include the name of the file in the material file (therefore in the undef_mat_prop)
  ! if we want more than one tomofile (Examples: 2 file with a differente resolution
  ! as in los angeles case we need to loop over mat_ext_mesh(1,ispec)...
  ! it is a possible solution )
  !  magnoni 1/12/09
  filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//trim(TOMO_FILENAME)
  open(unit=27,file=trim(filename),status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error reading tomography file')

  ! reads in model dimensions
  read(27,*) ORIG_X, ORIG_Y, ORIG_Z, END_X, END_Y, END_Z
  read(27,*) SPACING_X, SPACING_Y, SPACING_Z
  read(27,*) NX, NY, NZ
  read(27,*) VP_MIN, VP_MAX, VS_MIN, VS_MAX, RHO_MIN, RHO_MAX

  print *,'read NX',NX,NY,NZ,SPACING_X, SPACING_Y, SPACING_Z

  allocate(array_read_tmp(NX))
  nrecord = NX*NY*NZ

  ! allocates model records
  allocate(vp_tomography(1:nrecord), &
          vs_tomography(1:nrecord), &
          rho_tomography(1:nrecord), &
          z_tomography(1:nrecord),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  ! reads in record sections
!  do irecord = 1,nrecord
!    read(27,*) x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo

    ! stores record values
!    vp_tomography(irecord) = vp_tomo
!    vs_tomography(irecord) = vs_tomo
!    rho_tomography(irecord) = rho_tomo
!    z_tomography(irecord) = z_tomo
!  enddo
  do iz=1,NZ
    z_tmp=ORIG_Z+(iz-1)*SPACING_Z
    do iy=1,NY
      irecord0=(iz-1)*NX*NY+(iy-1)*NX+1
      read(27,*)array_read_tmp(:)
      vp_tomography(irecord0:irecord0+NX-1)=array_read_tmp(1:NX)
      z_tomography(irecord0:irecord0+NX-1)=z_tmp
    end do
  end do


  close(27)
  deallocate(array_read_tmp)

  ! user output
  if( myrank == 0 ) then
    write(IMAIN,*) '     tomography model: ',trim(TOMO_FILENAME)
  endif

  end subroutine read_model_tomography

  subroutine read_topography(myrank)
  use tomography

  implicit none

  integer :: myrank
  integer ::ilat,ilon
  character(len=256):: filename
  integer :: irecord,ier
  real(kind=CUSTOM_REAL) ::temp

  filename = MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES))//'SA_realtopo.dat'
  open(unit=27,file=trim(filename),status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error reading tomography file')

  ! reads in model dimensions
!  read(27,*) NLON, NLAT
!  read(27,*) ORIG_LON,ORIG_LAT, END_LON, END_LAT
!  read(27,*) SPACING_LON, SPACING_LAT
  NLON=3312        
  NLAT=2352 
  ORIG_LON=-85.7
  ORIG_LAT=-39.84 
  SPACING_LON=0.0083   
  SPACING_LAT=0.0083
  allocate(topography(NLAT,NLON))
!  do ilon=1,NLON
     do ilat=NLAT,1,-1
        read(27,*)topography(ilat,:)
     end do
!  end do


  end subroutine read_topography

  subroutine read_trench(myrank)
     use tomography

  implicit none

  integer :: myrank
  integer ::ilat,ilon
  character(len=256):: filename
  integer :: itrench,ier
  
  filename = MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES))//'trench.dat'
  open(unit=27,file=trim(filename),status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error reading tomography file')

  ! reads in model dimensions
!  read(27,*) NLON, NLAT
!  read(27,*) ORIG_LON,ORIG_LAT, END_LON, END_LAT
!  read(27,*) SPACING_LON, SPACING_LAT
  read(27,*)NTRENCH
  allocate(lat_trench(NTRENCH))
  allocate(lon_trench(NTRENCH))
!  do ilon=1,NLON
     do itrench=1,NTRENCH
        read(27,*)lon_trench(itrench),lat_trench(itrench)
     end do
!  end do
  
  end subroutine read_trench

!
!-------------------------------------------------------------------------------------------------
!
  subroutine  read_DSM_model (myrank)
   use tomography
!  use tomography, only:n_structure_zone,rmin_structure_zone,rmax_structure_zone,&
!               rho_structure_zone,vpv_structure_zone,vph_structure_zone,vsv_structure_zone,&
!               vsh_structure_zone,eta_structure_zone,qmu_structure_zone,qkappa_structure_zone,&
!               fluid_thiszone,DSM_FILENAME,CUSTOM_REAl
  implicit none

  integer ::myrank
!local parameters for reading file
  real(kind=CUSTOM_REAL),parameter ::kmtom=1000.0
  real(kind=CUSTOM_REAL),parameter ::gcmcubetokmmcube=1000.0
  real(kind=CUSTOM_REAL):: time_series_length,omega_imag
  integer:: n_frequency,ngrid_r,lmin,lmax
  real(kind=CUSTOM_REAL):: r_freesurf,r_ICB,r_CMB
  real(kind=CUSTOM_REAL)::source_r,source_depth,source_lat,source_lon
  integer ::source_type,nexp,save_velo
  real(kind=CUSTOM_REAL)::source_mt(3,3)
  real(kind=CUSTOM_REAL)::fr,ftheta,fphi
  character(len=80) ::DSM_file_name,tmp_file_name

  integer ::i
  

  DSM_file_name = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//&
                  trim(DSM_FILENAME)
! opening the temporary file
        open( unit=11, file=DSM_file_name, status='unknown' )
! reading the parameters
! ---- parameters for time series ---
        read(11,*) time_series_length,n_frequency
        read(11,*) omega_imag
! ---- parameters for numerical grids ----
        read(11,*) ngrid_r,lmin,lmax
! ---- parameters for structure ---
        read(11,*) n_structure_zone

!allocate arrays
        allocate(rmin_structure_zone(n_structure_zone))
        allocate(rmax_structure_zone(n_structure_zone))

        allocate(rho_structure_zone(4,n_structure_zone))
        allocate(vpv_structure_zone(4,n_structure_zone))
        allocate(vph_structure_zone(4,n_structure_zone))
        allocate(vsv_structure_zone(4,n_structure_zone))
        allocate(vsh_structure_zone(4,n_structure_zone))
        allocate(eta_structure_zone(4,n_structure_zone))
        allocate(qmu_structure_zone(n_structure_zone))
        allocate(qkappa_structure_zone(n_structure_zone))
        allocate(fluid_thiszone(n_structure_zone))
        r_CMB=0.d0
        r_ICB=0.d0
        do 130 i=1,n_structure_zone
          read(11,*) rmin_structure_zone(i),rmax_structure_zone(i),&
                     rho_structure_zone(1,i),rho_structure_zone(2,i),&
                     rho_structure_zone(3,i),rho_structure_zone(4,i)
          read(11,*) vpv_structure_zone(1,i),vpv_structure_zone(2,i),&
                     vpv_structure_zone(3,i),vpv_structure_zone(4,i)
          read(11,*) vph_structure_zone(1,i),vph_structure_zone(2,i),&
                     vph_structure_zone(3,i),vph_structure_zone(4,i)
          read(11,*) vsv_structure_zone(1,i),vsv_structure_zone(2,i),&
                     vsv_structure_zone(3,i),vsv_structure_zone(4,i)
          read(11,*) vsh_structure_zone(1,i),vsh_structure_zone(2,i),&
                     vsh_structure_zone(3,i),vsh_structure_zone(4,i)
          read(11,*) eta_structure_zone(1,i),eta_structure_zone(2,i),&
                     eta_structure_zone(3,i),eta_structure_zone(4,i),&
                     qmu_structure_zone(i),qkappa_structure_zone(i)

          rmin_structure_zone(i)=rmin_structure_zone(i)*kmtom
          rmax_structure_zone(i)=rmax_structure_zone(i)*kmtom
          rho_structure_zone(:,i)=rho_structure_zone(:,i)*gcmcubetokmmcube
          vpv_structure_zone(:,i)=vpv_structure_zone(:,i)*kmtom
          vph_structure_zone(:,i)=vph_structure_zone(:,i)*kmtom
          vsv_structure_zone(:,i)=vsv_structure_zone(:,i)*kmtom
          vsh_structure_zone(:,i)=vsh_structure_zone(:,i)*kmtom


          if ( ( ( vsv_structure_zone(1,i).eq.0.d0 ).and. &
                ( vsv_structure_zone(2,i).eq.0.d0 ).and. &
                ( vsv_structure_zone(3,i).eq.0.d0 ).and. &
                ( vsv_structure_zone(4,i).eq.0.d0 )      ).or.&
              ( ( vsh_structure_zone(1,i).eq.0.d0 ).and.&
                ( vsh_structure_zone(2,i).eq.0.d0 ).and.&
                ( vsh_structure_zone(3,i).eq.0.d0 ).and.&
                ( vsh_structure_zone(4,i).eq.0.d0 )    ) ) then
                fluid_thiszone(i)=1
          else
                fluid_thiszone(i)=0
          end if

          if(i.ge.2.and.i.lt.n_structure_zone) then
              if(abs(vsv_structure_zone(1,i)).lt.1.e-7.and.&
                abs(vsv_structure_zone(1,i-1)).gt.1.e-7) then
                   r_ICB=rmin_structure_zone(i)
              end if
             if(abs(vsv_structure_zone(1,i)).gt.1.e-7.and.&
                abs(vsv_structure_zone(1,i-1)).lt.1.e-7) then
                  r_CMB=rmin_structure_zone(i)
             end if

          end if
  130   continue
        r_freesurf=rmax_structure_zone(n_structure_zone)
        if(r_CMB.lt.1.e-7.or.r_ICB.lt.1.e-7) then
               print *,'r_CMB',r_CMB,r_ICB
               stop 'Error in finding r_CMB or r_ICB'
        end if

! ---- parameters for a source ---
        read(11,*) source_depth,source_lat,source_lon,source_type
! WENBO
        if(source_type.eq.1) then
           read(11,*) nexp,&
                  source_mt(1,1),source_mt(2,2),source_mt(3,3),&
                  source_mt(1,2),source_mt(1,3),source_mt(2,3)
        else if(source_type.eq.2) then
           read(11,*) nexp,fr,ftheta,fphi
        else
           stop 'Error of source_type'
        end if
! WENBO
        source_r = rmax_structure_zone(n_structure_zone) - source_depth
        source_mt(1,1) = source_mt(1,1) * ( 10.d0**(nexp-25) )
        source_mt(2,2) = source_mt(2,2) * ( 10.d0**(nexp-25) )
        source_mt(3,3) = source_mt(3,3) * ( 10.d0**(nexp-25) )
        source_mt(1,2) = source_mt(1,2) * ( 10.d0**(nexp-25) )
        source_mt(1,3) = source_mt(1,3) * ( 10.d0**(nexp-25) )
        source_mt(2,3) = source_mt(2,3) * ( 10.d0**(nexp-25) )
        source_mt(2,1) = source_mt(1,2)
        source_mt(3,1) = source_mt(1,3)
        source_mt(3,2) = source_mt(2,3)
        fr=fr*( 10.d0**(nexp-25) )
        ftheta=ftheta*( 10.d0**(nexp-25) )
        fphi=fphi*( 10.d0**(nexp-25) )

! --- parameters for stations ---
        read(11,*)tmp_file_name
        read(11,*)tmp_file_name
        read(11,*)tmp_file_name
        read(11,*)tmp_file_name

        read(11,*)save_velo
        close(11)

        print *,'read done',myrank
  end subroutine read_DSM_model


  subroutine model_tomography(flag_media,x_eval,y_eval,z_eval,r_middle, &
                             rho_final,vp_final,vs_final,myrank)

  use tomography

  implicit none

  double precision, intent(in) :: x_eval,y_eval,z_eval,r_middle
  real(kind=CUSTOM_REAL), intent(out) :: vp_final,vs_final,rho_final
  integer,intent(in) ::myrank,flag_media

  ! local parameters
  real(kind=CUSTOM_REAL) ::Radi,temp
  integer ::izone,izone1_Radi,izone2_Radi
  integer ::fluid_thispoint

  !sediment
  real(kind=CUSTOM_REAL) ::theta_surf,phi_surf
  integer ::ilat,ilon
  real(kind=CUSTOM_REAL) ::topo

!set up sediment close to trench
  real(kind=CUSTOM_REAL) ::dist_min,dist,azim,bazim
  real(kind=CUSTOM_REAL) ::dist_Ocbot,dist_Ocbot_square,dist_Ocbot_cube
  integer ::itrench,itrench_find

  Radi=dsqrt(x_eval**2+y_eval**2+z_eval**2)
  
!  if(ELLIPTICITY) then
!     dcost = dcos(ReceiverInfo(i)%theta)
!     p20 = 0.5d0*(3.0d0*dcost*dcost-1.0d0)
!     call spline_evaluation(rspl,espl,espl2,nspl,R_EARTH,ell)
!     Rearth  =  R_EARTH*(1.0d0-(2.0d0/3.0d0)*ell*p20)
!     lat= 180.0*PI*atan(tan(PI/2-theta)/0.99329534d0)
!  else  
!     Rearth  =  R_EARTH
!     lat= 180/PI*(PI/2-theta)
!  end if
!  if(abs(Radi-r_layer(2)).lt.100) then
!    if(myrank.eq.1) print *,'loca',Radi,r_middle

!print *,"fluid or not",flag_media
!solid media
  if(flag_media.eq.2) then
      fluid_thispoint=0
!fluid media
  else if(flag_media.eq.1) then
      fluid_thispoint=1
  else
      stop 'fluid_thispoint shoud be either 1(acoustic) or 2(elastic).'
  end if

  izone2_Radi=0
  izone1_Radi=1

  do izone=1,n_structure_zone
     if(rmin_structure_zone(izone).le.Radi) then
         izone1_Radi=izone
     end if
  end do

  do izone=n_structure_zone,1,-1
     if(abs(rmax_structure_zone(izone)-Radi).le.depth_tolerence) then
        izone1_Radi=izone
        if(izone.lt.n_structure_zone) then
          izone2_Radi=izone+1
!          print *,'izone2',izone2_Radi,Radi
          if(fluid_thispoint.ne.fluid_thiszone(izone2_Radi)) izone2_Radi=0
        end if
     end if
  end do
if(izone2_Radi.gt.n_structure_zone.and.fluid_thispoint.ne.fluid_thiszone(izone1_Radi))&
print *,'izone22',Radi,izone2_Radi,izone1_Radi,fluid_thispoint
  if(izone2_Radi.gt.n_structure_zone) izone2_Radi=0
!WENBO to be canceled
  if(fluid_thispoint.ne.fluid_thiszone(izone1_Radi)) &
           print *,"izone1",izone1_Radi,fluid_thiszone(izone1_Radi),fluid_thispoint
  if(fluid_thispoint.ne.fluid_thiszone(izone1_Radi)) izone1_Radi=0

!WENBO to be fixed and removed (canceled)
!Bathymetry
!Gofar fault
!  if(izone2_Radi.eq.0.and.izone1_Radi.eq.0.and.abs(Radi-6368000.0).lt.1900.0) then
!     if(fluid_thispoint.eq.1) then
!        izone2_Radi=11
!     else
!        izone1_Radi=10
!     end if
!     Radi=6368000.0
!  end if

!South America
  if(izone2_Radi.eq.0.and.izone1_Radi.eq.0.and.abs(Radi-6367122.0).lt.12000.0) then
     if(fluid_thispoint.eq.1) then
        izone2_Radi=11
     else
        izone1_Radi=10
     end if
     Radi=6367122.0
  end if


  if(izone2_Radi.eq.0.and.izone1_Radi.eq.0) then
       print *,'izone2',izone2_Radi,Radi,fluid_thispoint,rmax_structure_zone(1),&
                   r_middle,x_eval,y_eval,z_eval,180.0/3.14159*asin(z_eval/Radi),&
                   180.0/3.14159*atan(y_eval/(x_eval+1.e-9))
       stop 'Error, fluid_thispoint is not equal to fluid_izone found'
  else if(izone2_Radi.gt.0.and.izone1_Radi.gt.0) then
     if(r_middle.lt.Radi) then
       call cal_PREM_structure(Radi,rho_structure_zone(1,izone1_Radi),temp)
       rho_final=temp
       call cal_PREM_structure(Radi,vpv_structure_zone(1,izone1_Radi),temp)
       vp_final=temp
       call cal_PREM_structure(Radi,vsv_structure_zone(1,izone1_Radi),temp)
       vs_final=temp
     else 
       call cal_PREM_structure(Radi,rho_structure_zone(1,izone2_Radi),temp)
       rho_final=temp
       call cal_PREM_structure(Radi,vpv_structure_zone(1,izone2_Radi),temp)
       vp_final=temp
       call cal_PREM_structure(Radi,vsv_structure_zone(1,izone2_Radi),temp)
       vs_final=temp
     end if

!       if(myrank.eq.1) &
! for deep coupling check
!       print *,'structure1',6371000.0-Radi,6371000.0-r_middle,vp_final,vs_final,rho_final, &
!          izone1_Radi,izone2_Radi,rmin_structure_zone(izone1_Radi)

!To be removed????
!       if(Radi>6371000.0-9000.0+100.0.or.Radi<6371000.0-9000.0-100.0) &
!           stop 'Error discontinuity'

  else if(izone1_Radi.gt.0) then
       call cal_PREM_structure(Radi,rho_structure_zone(1,izone1_Radi),temp)
       rho_final=temp
       call cal_PREM_structure(Radi,vpv_structure_zone(1,izone1_Radi),temp)
       vp_final=temp
       call cal_PREM_structure(Radi,vsv_structure_zone(1,izone1_Radi),temp)
       vs_final=temp
!       if(Radi>6371000.0-9000.0+100.0.and.izone1_Radi.eq.9) stop 'Error izone9'
!       if(Radi<6371000.0-9000.0-100.0.and.izone1_Radi.eq.10) stop 'Error izone10'
 ! print *,'structure1',6371000.0-Radi,vp_final,vs_final,rho_final, &
 !         izone1_Radi,izone2_Radi,rmin_structure_zone(izone1_Radi)

  else if(izone2_Radi.gt.0) then
!       if(Radi>6371000.0-9000.0+100.0.and.izone2_Radi.eq.9) stop 'Error izone9'
!       if(Radi<6371000.0-9000.0-100.0.and.izone2_Radi.eq.10) stop 'Error izone10'
       call cal_PREM_structure(Radi,rho_structure_zone(1,izone2_Radi),temp)
       rho_final=temp
       call cal_PREM_structure(Radi,vpv_structure_zone(1,izone2_Radi),temp)
       vp_final=temp
       call cal_PREM_structure(Radi,vsv_structure_zone(1,izone2_Radi),temp)
       vs_final=temp
!        print *,'izone2_error',izone2_Radi,Radi,vp_final,fluid_thispoint
!        stop 'other errors_zone2'
   else
        stop 'other errors'
  end if

!South America - sediment

!Radi might be changed, so recalculate it.
              Radi=dsqrt(x_eval**2+y_eval**2+z_eval**2)
              theta_surf=dacos(z_eval/Radi)
!            sin(phi_surf)=ystore(iglob)/r_surf/dsin(theta_surf)
              phi_surf=dacos(dble(x_eval)/(dble(Radi)*dsin(theta_surf))*(1-1.e-7))
              if(y_eval<-1.e-14) then
                    phi_surf=-phi_surf
              end if
              theta_surf=theta_surf*180.0/3.1415926
              phi_surf=phi_surf*180.0/3.1415926
              ilat=(90.0-theta_surf-ORIG_LAT)/SPACING_LAT
              ilon=(phi_surf-ORIG_LON)/SPACING_LON
              topo=topography(ilat,ilon)

              dist_min=180.0
              do itrench=1,NTRENCH
                 call  distazbaz(90.0-theta_surf,phi_surf,lat_trench(itrench),&
                       lon_trench(itrench),dist,azim,bazim)
                       if(dist<dist_min) then
                          dist_min=dist
                          itrench_find=itrench
                       end if
              end do
!              if(myrank.eq.55) print *,'dist=',dist_min,90.0-theta_surf,phi_surf,lat_trench(itrench_find),&
!                lon_trench(itrench_find)
!sediment layer. Vs, Vp and rho relationship is refered to Thomas M. Brocher,
!BSSA(2010)
              if(Radi-R_EARTH_SURF-topo.gt.-3000.0.and.&
                 topo.lt.0.and.fluid_thispoint.eq.0&
                 .and.phi_surf.ge.lon_trench(itrench_find)) then
                   vs_final=1000.0-(Radi-R_EARTH_SURF-topo)/3000.0*2000.0
                   if(vs_final.lt.1000.0) vs_final=1000.0
                   vs_final=vs_final/1000.0
                   vp_final=1.16*vs_final+1.36
                   if(vp_final.gt.4.0) then
                     vp_final=4.0
                     vs_final=(vp_final-1.36)/1.16
                   end if

                   rho_final=1.6612*vp_final-0.4721*vp_final**2+0.0671*vp_final**3&
                             -0.0043*vp_final**4+0.000106*vp_final**5
                   vs_final=vs_final*1000.0
                   vp_final=vp_final*1000.0
                   rho_final=rho_final*1000.0
              end if
!Marine  sediment
!Edwin L. Hmilton, JASA, 1979 and JSP, 1976
              if(Radi-R_EARTH_SURF-topo.gt.-2300.0.and.&
                 topo.lt.0.and.fluid_thispoint.eq.0&
                 .and.dist_min.le.0.15) then
                   dist_Ocbot=-(Radi-R_EARTH_SURF-topo)/1000.0
                   dist_Ocbot_square=dist_Ocbot*dist_Ocbot
                   dist_Ocbot_cube=dist_Ocbot_square*dist_Ocbot
                   vp_final=1.511+dist_Ocbot*1.304-0.741*dist_Ocbot_square+0.257*dist_Ocbot_cube
                   vs_final=0.78*vp_final-0.962
                   rho_final=1.53+1.395*dist_Ocbot-0.617*dist_Ocbot_square
                   if(vp_final.gt.4.0) then
                     vp_final=4.0
                     vs_final=0.78*vp_final-0.962
                   end if

                   vs_final=vs_final*1000.0
                   vp_final=vp_final*1000.0
                   rho_final=rho_final*1000.0

                   if(vs_final.lt.1000.0) vs_final=1000.0

              end if


  end subroutine model_tomography



  subroutine model_tomography1(x_eval,y_eval,z_eval,r_middle, &
                             rho_final,vp_final,vs_final,myrank)

  use tomography

  implicit none

  double precision, intent(in) :: x_eval,y_eval,z_eval,r_middle
  real(kind=CUSTOM_REAL), intent(out) :: vp_final,vs_final,rho_final
  integer ::myrank

  ! local parameters

  double precision ::rho_layer(4,3),vp_layer(4,3),vs_layer(4,3)
  double precision ::r_layer(3)
  double precision ::rearth,Radi,rend,temp
  integer ::ilayer

  rearth=6371000.d0
  r_layer(1)=6371000.d0
  rho_layer(1,1)=2900.d0
  rho_layer(2,1)=0.0
  rho_layer(3,1)=0.d0
  rho_layer(4,1)=0.d0
  vp_layer(1,1)=6500.0
  vp_layer(2,1)=0.0
  vp_layer(3,1)=0.d0
  vp_layer(4,1)=0.d0
  vs_layer(1,1)=3700.0
  vs_layer(2,1)=0.d0
  vs_layer(3,1)=0.d0
  vs_layer(4,1)=0.d0

  r_layer(2)=6362000.d0
  rho_layer(1,2)=2691.d0
  rho_layer(2,2)=692.4
  rho_layer(3,2)=0.d0
  rho_layer(4,2)=0.d0
  vp_layer(1,2)=8785.4
  vp_layer(2,2)=-749.5
  vp_layer(3,2)=0.d0
  vp_layer(4,2)=0.d0
  vs_layer(1,2)=4750.0
  vs_layer(2,2)=0.d0
  vs_layer(3,2)=0.d0
  vs_layer(4,2)=0.d0

  r_layer(3)=6271000.d0
  rho_layer(1,3)=2691.d0
  rho_layer(2,3)=692.4
  rho_layer(3,3)=0.d0
  rho_layer(4,3)=0.d0
  vp_layer(1,3)=25413.9
  vp_layer(2,3)=-17697.2
  vp_layer(3,3)=0.d0
  vp_layer(4,3)=0.d0
  vs_layer(1,3)=-12811361.4842722
  vs_layer(2,3)=40496077.735693
  vs_layer(3,3)=-42638043.8962742
  vs_layer(4,3)=14959075.7253793

  rend=6161000.d0


  !WENBO
  Radi=dsqrt(x_eval**2+y_eval**2+z_eval**2)
  if(Radi.le.r_layer(1)*(1.d0+1.e-7).and.Radi.gt.r_layer(2)) then
    ilayer=1
  else if(Radi.le.r_layer(2).and.Radi.gt.r_layer(3)) then
    ilayer=2
  else if(Radi.le.r_layer(3).and.Radi.ge.rend) then
    ilayer=3
  else
    print *,"Radi",Radi
    stop 'Error in setting vp,vs or rho'
  end if

  if(r_middle.gt.r_layer(2)) then
        ilayer=1
  else
        ilayer=2
  end if
!for deep coupling check
!  if(ilayer.ne.1) then
!     print *,'radius',ilayer,r_middle,radi
!  end if

  call cal_PREM_structure(Radi,rho_layer(1,ilayer),temp)
  rho_final=temp
  call cal_PREM_structure(Radi,vp_layer(1,ilayer),temp)
  vp_final=temp
  call cal_PREM_structure(Radi,vs_layer(1,ilayer),temp)
  vs_final=temp
if(abs(Radi-r_layer(2)).lt.100) then
  ilayer=1
  call cal_PREM_structure(Radi,rho_layer(1,ilayer),temp)
  rho_final=temp/2.0
  call cal_PREM_structure(Radi,vp_layer(1,ilayer),temp)
  vp_final=temp/2.0
  call cal_PREM_structure(Radi,vs_layer(1,ilayer),temp)
  vs_final=temp/2.0
  ilayer=2
  call cal_PREM_structure(Radi,rho_layer(1,ilayer),temp)
  rho_final=rho_final+temp/2.0
  call cal_PREM_structure(Radi,vp_layer(1,ilayer),temp)
  vp_final=vp_final+temp/2.0
  call cal_PREM_structure(Radi,vs_layer(1,ilayer),temp)
  vs_final=vs_final+temp/2.0

end if
!  print *,'structure',6371000.0-Radi,vp_final,vs_final,rho_final
  if(vs_final.lt.1.0) stop 'vs error'


  end subroutine model_tomography1


 subroutine model_ningxia(x_eval,y_eval,z_eval,r_middle, &
                             rho_final,vp_final,vs_final,qkappa_atten,qmu_atten,imaterial_id)

  use tomography

  implicit none

  double precision, intent(in) :: x_eval,y_eval,z_eval,r_middle
  real(kind=CUSTOM_REAL), intent(out) :: vp_final,vs_final,rho_final
  real(kind=CUSTOM_REAL), intent(out) :: qkappa_atten,qmu_atten
  integer, intent(in) :: imaterial_id

  integer ::myrank

  ! local parameters

  double precision ::rho_layer(4,3),vp_layer(4,3),vs_layer(4,3)
  double precision ::r_layer(3)
  double precision ::rearth,Radi,rend,temp
  integer ::ilayer

  qmu_atten = 80.0
  qkappa_atten=9999.

  if(z_eval.gt.-30000.0) then
     vp_final=6800.0
     vs_final=3900.0
     rho_final=2900.0
  else if(z_eval.gt.-101000.0) then
     vp_final=8110.0
     vs_final=4490.0
     rho_final=3380.0
  else 
     stop 'Error tomography_ningxia'
  end if



  call get_perturbation(x_eval,y_eval,z_eval, &
                             rho_final,vp_final,vs_final)
  !print
  !*,'structure',6371000.0-Radi,vp_final,vs_final,rho_final,x_eval,y_eval,z_eval
!  if(vs_final.lt.1.0) stop 'vs error'


  end subroutine model_ningxia





  subroutine cal_PREM_structure(r,param,final_value)
   use constants
   implicit none
   !double precision, intent(in)::r
   !double precision, intent(in)::param(4)
   !double precision, intent(out) ::final_value
   !double precision::a,rearth_prem

   real(kind=CUSTOM_REAL), intent(in)::r
   real(kind=CUSTOM_REAL), intent(in)::param(4)
   real(kind=CUSTOM_REAL), intent(out) ::final_value
   real(kind=CUSTOM_REAL)::a,rearth_prem
   rearth_prem=6371000.0
   a=r/rearth_prem
   final_value=param(1) &
               + param(2)*a &
               + param(3)*a*a &
               + param(4)*a*a*a  
   if(final_value.lt.1.d0) then
   !  print *,'par',param(:),a
   end if

  end subroutine cal_PREM_structure


!***********************find the perturbation***********************
  subroutine get_perturbation(x_eval,y_eval,z_eval, &
                             rho_final,vp_final,vs_final)
  use tomography

  implicit none

  !integer, intent(in) :: NX,NY,NZ
  !real(kind=CUSTOM_REAL), dimension(1:NX*NY*NZ), intent(in) ::
  !vp_tomography,vs_tomography,rho_tomography,z_tomography
  !double precision, intent(in) ::
  !ORIG_X,ORIG_Y,ORIG_Z,SPACING_X,SPACING_Y,SPACING_Z
  !double precision, intent(in) :: VP_MIN,VS_MIN,RHO_MIN,VP_MAX,VS_MAX,RHO_MAX

  double precision, intent(in) :: x_eval,y_eval,z_eval
  real(kind=CUSTOM_REAL), intent(out) :: vp_final,vs_final,rho_final

  ! local parameters
  integer :: ix,iy,iz
  integer :: p0,p1,p2,p3,p4,p5,p6,p7

  double precision ::x_eval_proj,y_eval_proj,z_eval_proj
  double precision :: spac_x,spac_y,spac_z
  double precision :: gamma_interp_x,gamma_interp_y
  double precision :: gamma_interp_z1,gamma_interp_z2,gamma_interp_z3, &
    gamma_interp_z4,gamma_interp_z5,gamma_interp_z6,gamma_interp_z7,gamma_interp_z8
  real(kind=CUSTOM_REAL) :: vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8, &
    vs1,vs2,vs3,vs4,vs5,vs6,vs7,vs8,rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8
  real(kind=CUSTOM_REAL) :: perturbation
  logical ::abs_VpVsRho


  double precision ::Radi,theta,Phi,depth,lat,lon,Rearth

  !WENBO
  Radi=dsqrt(x_eval**2+y_eval**2+z_eval**2)
!  theta=dacos(z_eval/Radi)
  if(y_eval<0.0) then
!       Phi=2*PI-dacos(x_eval*(1.d0-1.e-10)/(Radi*dsin(Theta)))
  else
!       Phi=dacos(x_eval*(1.d0-1.e-10)/(Radi*dsin(Theta)))
  end if
!  x_eval_proj=Radi*dcos(Phi)
!  y_eval_proj=Radi*dsin(Phi)
!  z_eval_proj=0.0
  x_eval_proj=x_eval
  y_eval_proj=y_eval
  z_eval_proj=z_eval

  !lon=Phi/PI*180.0

!  if(ELLIPTICITY) then
!     dcost = dcos(ReceiverInfo(i)%theta)
!     p20 = 0.5d0*(3.0d0*dcost*dcost-1.0d0)
!     call spline_evaluation(rspl,espl,espl2,nspl,R_EARTH,ell)
!     Rearth  =  R_EARTH*(1.0d0-(2.0d0/3.0d0)*ell*p20)
!     lat= 180.0*PI*atan(tan(PI/2-theta)/0.99329534d0)
!  else  
!     Rearth  =  R_EARTH
!     lat= 180/PI*(PI/2-theta)
!  end if

  Rearth=R_TOP_BOUND
  !lat= 180.0/PI*(PI/2-theta)

  depth=Rearth-Radi
  ! determine spacing and cell for linear interpolation
  
  !spac_x = (lon - ORIG_X) / SPACING_X
  !spac_y = (lat - ORIG_Y) / SPACING_Y
  !spac_z = (depth - ORIG_Z) / SPACING_Z
  spac_x = (x_eval_proj - ORIG_X) / SPACING_X
  spac_y = (y_eval_proj - ORIG_Y) / SPACING_Y
  spac_z = (z_eval_proj - ORIG_Z) / SPACING_Z

  !WENBO
   !print *,'model',lon,lat,depth,ORIG_X,ORIG_Y,ORIG_Z
  ix = int(spac_x)
  iy = int(spac_y)
  iz = int(spac_z)

!for deep coupling check
!  if(ix.lt.1.or.ix.gt.NX.or.iy.lt.1.or.iy.gt.NY.or.iz.lt.1.or.iz.gt.NZ) &
!         print *,'ixiyiz',ix,iy,iz


  gamma_interp_x = spac_x - dble(ix)
  gamma_interp_y = spac_y - dble(iy)

  ! suppress edge effects for points outside of the model SPOSTARE DOPO
  if(ix < 0) then
    ix = 0
    gamma_interp_x = 0.d0
  endif
  if(ix > NX-2) then
    ix = NX-2
    gamma_interp_x = 1.d0
  endif

  if(iy < 0) then
    iy = 0
    gamma_interp_y = 0.d0
  endif
  if(iy > NY-2) then
    iy = NY-2
    gamma_interp_y = 1.d0
  endif

  if(iz < 0) then
     iz = 0
  !   gamma_interp_z = 0.d0
  endif
  if(iz > NZ-2) then
     iz = NZ-2
  !  gamma_interp_z = 1.d0
  endif


  ! define 8 corners of interpolation element
  p0 = ix+iy*NX+iz*(NX*NY)
  p1 = (ix+1)+iy*NX+iz*(NX*NY)
  p2 = (ix+1)+(iy+1)*NX+iz*(NX*NY)
  p3 = ix+(iy+1)*NX+iz*(NX*NY)
  p4 = ix+iy*NX+(iz+1)*(NX*NY)
  p5 = (ix+1)+iy*NX+(iz+1)*(NX*NY)
  p6 = (ix+1)+(iy+1)*NX+(iz+1)*(NX*NY)
  p7 = ix+(iy+1)*NX+(iz+1)*(NX*NY)

  gamma_interp_z1=0.d0
  gamma_interp_z2=0.d0
  gamma_interp_z3=0.d0
  gamma_interp_z4=0.d0

 

  if(z_tomography(p4+1) == z_tomography(p0+1)) then
          gamma_interp_z1 = 1.d0
  else
     if(abs(z_tomography(p4+1)-z_tomography(p0+1)).lt.1.e-5) &
      print *,'floating1',z_tomography(p4+1),z_tomography(p0+1),p4,p0
          gamma_interp_z1 = (z_eval-z_tomography(p0+1))/(z_tomography(p4+1)-z_tomography(p0+1))
  endif
  if(gamma_interp_z1 > 1.d0) then
          gamma_interp_z1 = 1.d0
  endif
  if(gamma_interp_z1 < 0.d0) then
          gamma_interp_z1 = 0.d0
  endif


  if(z_tomography(p5+1) == z_tomography(p1+1)) then
          gamma_interp_z2 = 1.d0
  else
     if(abs(z_tomography(p5+1)-z_tomography(p1+1)).lt.1.e-5) &
      print *,'floating2',z_tomography(p5+1),z_tomography(p1+1),p5,p1
          gamma_interp_z2 = (z_eval-z_tomography(p1+1))/(z_tomography(p5+1)-z_tomography(p1+1))
  endif
  if(gamma_interp_z2 > 1.d0) then
          gamma_interp_z2 = 1.d0
  endif
  if(gamma_interp_z2 < 0.d0) then
          gamma_interp_z2 = 0.d0
  endif


  if(z_tomography(p6+1) == z_tomography(p2+1)) then
          gamma_interp_z3 = 1.d0
  else
     if(abs(z_tomography(p6+1)-z_tomography(p2+1)).lt.1.e-5) &
      print *,'floating3',z_tomography(p6+1),z_tomography(p2+1),p6,p2
          gamma_interp_z3 = (z_eval-z_tomography(p2+1))/(z_tomography(p6+1)-z_tomography(p2+1))
  endif
  if(gamma_interp_z3 > 1.d0) then
          gamma_interp_z3 = 1.d0
  endif
  if(gamma_interp_z3 < 0.d0) then
          gamma_interp_z3 = 0.d0
  endif


  if(z_tomography(p7+1) == z_tomography(p3+1)) then
          gamma_interp_z4 = 1.d0
  else
     if(abs(z_tomography(p7+1)-z_tomography(p3+1)).lt.1.e-5) &
      print *,'floating4',z_tomography(p7+1),z_tomography(p3+1),p7,p3
          gamma_interp_z4 = (z_eval-z_tomography(p3+1))/(z_tomography(p7+1)-z_tomography(p3+1))
  endif
  if(gamma_interp_z4 > 1.d0) then
          gamma_interp_z4 = 1.d0
  endif
  if(gamma_interp_z4 < 0.d0) then
          gamma_interp_z4 = 0.d0
  endif

  gamma_interp_z5 = 1. - gamma_interp_z1
  gamma_interp_z6 = 1. - gamma_interp_z2
  gamma_interp_z7 = 1. - gamma_interp_z3
  gamma_interp_z8 = 1. - gamma_interp_z4

  vp1 = vp_tomography(p0+1)
  vp2 = vp_tomography(p1+1)
  vp3 = vp_tomography(p2+1)
  vp4 = vp_tomography(p3+1)
  vp5 = vp_tomography(p4+1)
  vp6 = vp_tomography(p5+1)
  vp7 = vp_tomography(p6+1)
  vp8 = vp_tomography(p7+1)

  vs1 = vs_tomography(p0+1)
  vs2 = vs_tomography(p1+1)
  vs3 = vs_tomography(p2+1)
  vs4 = vs_tomography(p3+1)
  vs5 = vs_tomography(p4+1)
  vs6 = vs_tomography(p5+1)
  vs7 = vs_tomography(p6+1)
  vs8 = vs_tomography(p7+1)

  rho1 = rho_tomography(p0+1)
  rho2 = rho_tomography(p1+1)
  rho3 = rho_tomography(p2+1)
  rho4 = rho_tomography(p3+1)
  rho5 = rho_tomography(p4+1)
  rho6 = rho_tomography(p5+1)
  rho7 = rho_tomography(p6+1)
  rho8 = rho_tomography(p7+1)

  ! use trilinear interpolation in cell to define Vp Vs and rho
  abs_VpVsRho=.false.
  if(abs_VpVsRho) then
   vp_final = &
     vp1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + &
     vp2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + &
     vp3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + &
     vp4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + &
     vp5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + &
     vp6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + &
     vp7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + &
     vp8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4

   vs_final = &
     vs1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + &
     vs2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + &
     vs3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + &
     vs4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + &
     vs5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + &
     vs6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + &
     vs7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + &
     vs8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4

   rho_final = &
     rho1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + &
     rho2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + &
     rho3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + &
     rho4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + &
     rho5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + &
     rho6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + &
     rho7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + &
     rho8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4
!WENBO
!     vs_final=vp_final/1.732
!     rho_final = (1.6612 * (vp_final / 1000.0)     &
!                      -0.4720 * (vp_final / 1000.0)**2  &
!                      +0.0671 * (vp_final / 1000.0)**3  &
!                      -0.0043 * (vp_final / 1000.0)**4  &
!                      +0.000106*(vp_final / 1000.0)**5 )*1000.0


  ! impose minimum and maximum velocity and density if needed
     if(vp_final < VP_MIN) vp_final = VP_MIN
     if(vs_final < VS_MIN) vs_final = VS_MIN
     if(rho_final < RHO_MIN) rho_final = RHO_MIN
     if(vp_final > VP_MAX) vp_final = VP_MAX
     if(vs_final > VS_MAX) vs_final = VS_MAX
     if(rho_final > RHO_MAX) rho_final = RHO_MAX
  
!WENBO
     if(ix.eq.(NX-2).or.ix.eq.0.or.iy.eq.NY-2.or.&
        iy.eq.0.or.iz.eq.NZ-2.or.iz.eq.0) then
         vp_final=vp1;vs_final=vs1;rho_final=rho1
     end if
  else
     perturbation= &
       vp1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + &
       vp2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + &
       vp3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + &
       vp4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + &
       vp5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + &
       vp6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + &
       vp7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + &
       vp8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4
     if(vs_final.gt.900.0) then
       vp_final=vp_final*(1+0.5*perturbation)
       vs_final=vs_final*(1+1.0*perturbation)
!       rho_final=rho_final*(1+2.0*perturbation)

!       if(z_eval.gt.-50.0.and.z_eval.lt.50.0) print *,"vpvsrho",x_eval,y_eval,perturbation 
!       if(x_eval.gt.5550050.0.and.x_eval.lt.5550100.0)  &
!           write(*,"(A,3X,3E14.7)") "vpvsrho1",z_eval,y_eval,perturbation

       if(perturbation.gt.0.2.or.perturbation.lt.-0.2) stop 'Error perturbation'
     end if

!    vp_final=vp1
!    vs_final=vp1/1.732
!    rho_final = (1.6612 * (vp_final / 1000.0)     &
!                  -0.4720 * (vp_final / 1000.0)**2  &
!                  +0.0671 * (vp_final / 1000.0)**3  &
!                  -0.0043 * (vp_final / 1000.0)**4  &
!                  +0.000106*(vp_final / 1000.0)**5 )*1000.0

   !print *,lat,lon,depth,vp_final,vs_final,vs1,vs8,rho_final,ix,iy,iz,NX,NY,NZ
  end if
  !if(depth<30000) then
  !    vp_final=5800.0;vs_final=3200.0;rho_final=2600.0
  !else 
  !    vp_final=8100.0;vs_final=4484.0;rho_final=3379.0
  !end if

  end subroutine get_perturbation


  subroutine deallocate_tomography_files()

    use tomography

    implicit none

    ! deallocates models dimensions
    !deallocate(ORIG_X,ORIG_Y,ORIG_Z)
    !deallocate(SPACING_X,SPACING_Y,SPACING_Z)

    ! deallocates models parameter records
    deallocate(vp_tomography)
    deallocate(vs_tomography)
    deallocate(rho_tomography)
    deallocate(z_tomography)

    ! deallocates models entries
    !deallocate(NX,NY,NZ)
    !deallocate(nrecord)

    ! deallocates models min/max statistics
    !deallocate(VP_MIN,VS_MIN,RHO_MIN)
    !deallocate(VP_MAX,VS_MAX,RHO_MAX)

    ! q values
!    if (any(tomo_has_q_values)) then
!      deallocate(qp_tomography)
!      deallocate(qs_tomography)
!    endif
!    deallocate(tomo_has_q_values)

  end subroutine deallocate_tomography_files

