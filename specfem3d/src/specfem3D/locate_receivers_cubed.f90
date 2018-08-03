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

!----
!---- locate_receivers finds the correct position of the receivers
!----
  subroutine locate_receivers_cubed(ibool,myrank,NSPEC_AB,NGLOB_AB,NGNOD, &
                 xstore,ystore,zstore,xigll,yigll,zigll,rec_filename, &
                 nrec,islice_selected_rec,ispec_selected_rec, &
                 xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
                 NPROC,theta_source,phi_source,CUBED_SPHERE_PROJECTION, &
                 iglob_is_surface_external_mesh,ispec_is_surface_external_mesh, &
                 num_free_surface_faces,free_surface_ispec,free_surface_ijk)
  use constants
  implicit none

!  include "constants.h"


  integer NPROC
  integer NGNOD
  logical CUBED_SPHERE_PROJECTION

  integer nrec,myrank,nrec_found

  integer NSPEC_AB,NGLOB_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

! for surface locating and normal computing with external mesh
!  integer :: pt0_ix,pt0_iy,pt0_iz,pt1_ix,pt1_iy,pt1_iz,pt2_ix,pt2_iy,pt2_iz
  integer :: num_free_surface_faces
  double precision r_surf,theta_surf,phi_surf
!  real(kind=CUSTOM_REAL), dimension(3) :: u_vector,v_vector,w_vector
  logical, dimension(NGLOB_AB) :: iglob_is_surface_external_mesh
  logical, dimension(NSPEC_AB) :: ispec_is_surface_external_mesh
  integer, dimension(num_free_surface_faces) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces) :: free_surface_ijk

  integer, allocatable, dimension(:) :: ix_initial_guess,iy_initial_guess,iz_initial_guess

  integer iprocloop
  integer iorientation
  double precision stazi,stdip
  integer ios

  double precision,dimension(1) :: altitude_rec,distmin_ele
  double precision,dimension(4) :: elevation_node,dist_node
  double precision,dimension(NPROC) :: distmin_ele_all,elevation_all
  double precision, allocatable, dimension(:) :: r_target
  double precision, allocatable, dimension(:) :: x_target,y_target,z_target
  double precision, allocatable, dimension(:) :: epidist
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision, allocatable, dimension(:,:) :: x_found_all,y_found_all,z_found_all

  integer irec
  integer i,j,k,ispec,iglob,iface,inode,imin,imax,jmin,jmax,kmin,kmax,igll,jgll,kgll
  integer iselected,jselected,iface_selected,iadjust,jadjust
  integer iproc(1)

  double precision n(3)
  double precision thetan,phin
  double precision sint,cost,sinp,cosp
!  double precision r0,p20
  double precision theta,phi
  double precision theta_source,phi_source
  double precision dist
  double precision xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

! input receiver file name
  character(len=256) rec_filename

! topology of the control points of the surface element
  integer iax,iay,iaz
  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! coordinates of the control points of the surface element
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer iter_loop,ispec_iterate

  integer ia
  double precision x,y,z
  double precision xix,xiy,xiz
  double precision etax,etay,etaz
  double precision gammax,gammay,gammaz

! timer MPI
  double precision, external :: wtime
  double precision time_start,tCPU

! use dynamic allocation
  double precision, dimension(:), allocatable :: final_distance
  double precision, dimension(:,:), allocatable :: final_distance_all
  double precision distmin,final_distance_max

! receiver information
! timing information for the stations
! station information for writing the seismograms

  integer :: iglob_selected
  integer, dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(3,3,nrec) :: nu
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name


  integer,allocatable, dimension(:) :: islice_selected_rec_found,ispec_selected_rec_found
  double precision,allocatable, dimension(:) :: xi_receiver_found,eta_receiver_found,gamma_receiver_found
  double precision,allocatable, dimension(:,:,:) :: nu_found
  character(len=MAX_LENGTH_STATION_NAME),allocatable, dimension(:) :: station_name_found
  character(len=MAX_LENGTH_NETWORK_NAME),allocatable, dimension(:) :: network_name_found
  double precision,allocatable,dimension(:) :: stlat_found,stlon_found,stele_found,stbur_found,epidist_found
!  character(len=150) STATIONS


  integer, allocatable, dimension(:,:) :: ispec_selected_rec_all
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur,elevation
  double precision, allocatable, dimension(:,:) :: xi_receiver_all,eta_receiver_all,gamma_receiver_all
  double precision, allocatable, dimension(:,:,:,:) :: nu_all


!  character(len=256) OUTPUT_FILES
  logical ::logical_tmp
!The below line is used to avoid possible error reports during compling
!the code. CUBED_SPHERE_PROJECTION is not used now, that occationally causes error 
!reports. But it might be usful in the future, so we keep them here.
  logical_tmp=CUBED_SPHERE_PROJECTION

  if(DEBUG_COUPLING) print *,'recv_file',rec_filename
  if(DEBUG_COUPLING) print *,'start',myrank
  ! get MPI starting time
  time_start = wtime()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
  endif

  ! define topology of the control element
  call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)
  !if(DEBUG_COUPLING) print *,'hex',myrank,rec_filename
  call synchronize_all()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '*****************************************************************'
    write(IMAIN,'(1x,a,a,a)') 'reading receiver information from ', trim(rec_filename), ' file'
    write(IMAIN,*) '*****************************************************************'
  endif
  if(DEBUG_COUPLING) print *,'readdata',myrank
  call synchronize_all()

  ! get number of stations from receiver file
  open(unit=1,file=trim(rec_filename),status='old',action='read',iostat=ios)
  if (ios /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))
  ! allocate memory for arrays using number of stations
  allocate(stlat(nrec))
  allocate(stlon(nrec))
  allocate(stele(nrec))
  allocate(stbur(nrec))
  allocate(epidist(nrec))
  allocate(elevation(nrec))

  allocate(ix_initial_guess(nrec))
  allocate(iy_initial_guess(nrec))
  allocate(iz_initial_guess(nrec))
  allocate(r_target(nrec))
  allocate(x_target(nrec))
  allocate(y_target(nrec))
  allocate(z_target(nrec))
  allocate(x_found(nrec))
  allocate(y_found(nrec))
  allocate(z_found(nrec))
  allocate(final_distance(nrec))

  allocate(islice_selected_rec_found(nrec))
  allocate(ispec_selected_rec_found(nrec))
  allocate(xi_receiver_found(nrec))
  allocate(eta_receiver_found(nrec))
  allocate(gamma_receiver_found(nrec))
  allocate(nu_found(3,3,nrec))
  allocate(station_name_found(nrec))
  allocate(network_name_found(nrec))
  allocate(stlat_found(nrec))
  allocate(stlon_found(nrec))
  allocate(stele_found(nrec))
  allocate(stbur_found(nrec))
  allocate(epidist_found(nrec))

  allocate(ispec_selected_rec_all(nrec,0:NPROC-1))
  allocate(xi_receiver_all(nrec,0:NPROC-1))
  allocate(eta_receiver_all(nrec,0:NPROC-1))
  allocate(gamma_receiver_all(nrec,0:NPROC-1))
  allocate(x_found_all(nrec,0:NPROC-1))
  allocate(y_found_all(nrec,0:NPROC-1))
  allocate(z_found_all(nrec,0:NPROC-1))
  allocate(final_distance_all(nrec,0:NPROC-1))
  allocate(nu_all(3,3,nrec,0:NPROC-1))

  ! loop on all the stations
  do irec=1,nrec

    read(1,*,iostat=ios) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
    if (ios /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))
    if(DEBUG_COUPLING) print *,'myrank',myrank,station_name(irec),network_name(irec),stlat(irec),stlon(irec)

    theta = PI/2.0d0 - stlat(irec)*PI/180.0d0
    !  if(.not.ELLIPTICITY) then
    !    theta = PI/2.0d0 - stlat(isource)*PI/180.0d0
    !  else
    !    theta = PI/2.0d0 - atan(0.99329534d0*dtan(stlat(isource)*PI/180.0d0))
    !  endif
    phi = stlon(irec)*PI/180.0d0
    if(myrank == 0.and.DEBUG_COUPLING) &
        print *,'myrank',myrank,'thetaphi',theta,phi,cos(theta)*cos(theta_source) + &
              sin(theta)*sin(theta_source)*cos(phi-phi_source)
    ! compute epicentral distance
    epidist(irec) = acos((cos(theta)*cos(theta_source) + &
              sin(theta)*sin(theta_source)*cos(phi-phi_source))*(1.0-TINYVAL))*180.0d0/PI
    if(myrank == 0.and.DEBUG_COUPLING) &
        print *,'myrank',myrank,'epidist(irec)',epidist(irec)



    ! print some information about stations
    if(myrank == 0) &
        write(IMAIN,*) 'Station #',irec,': ',station_name(irec)(1:len_trim(station_name(irec))), &
                       '.',network_name(irec)(1:len_trim(network_name(irec))), &
                       '    horizontal distance:  ',sngl(epidist(irec)),' deg'

    ! get approximate topography elevation at source long/lat coordinates
    !   set distance to huge initial value
    distmin = HUGEVAL
    if(num_free_surface_faces > 0) then
      iglob_selected = 1
      ! loop only on points inside the element
      ! exclude edges to ensure this point is not shared with other elements
      imin = 2
      imax = NGLLX - 1
      jmin = 2
      jmax = NGLLY - 1
      iselected = 0
      jselected = 0
      iface_selected = 0
      do iface=1,num_free_surface_faces
        do j=jmin,jmax
          do i=imin,imax

            ispec = free_surface_ispec(iface)
            igll = free_surface_ijk(1,(j-1)*NGLLY+i,iface)
            jgll = free_surface_ijk(2,(j-1)*NGLLY+i,iface)
            kgll = free_surface_ijk(3,(j-1)*NGLLY+i,iface)
            iglob = ibool(igll,jgll,kgll,ispec)
            r_surf=dsqrt(dble(xstore(iglob)**2)+dble(ystore(iglob)**2)+dble(zstore(iglob)**2))
            theta_surf=dacos(zstore(iglob)/r_surf)
!           sin(phi_surf)=ystore(iglob)/r_surf/dsin(theta_surf)
            phi_surf=dacos(dble(xstore(iglob))/(dble(r_surf)*dsin(theta_surf))*(1-1.0e-7))
            if(ystore(iglob)<-1.e-14) then
               phi_surf=2*PI-phi_surf
            end if
            ! keep this point if it is closer to the receiver
!            if(ELLIPTICITY) then
!                dcost = dcos(theta_surf)
!                p20 = 0.5d0*(3.0d0*dcost*dcost-1.0d0)
!                call spline_evaluation(rspl,espl,espl2,nspl,R_EARTH,ell)
!                r_ellip = R_EARTH*(1.0d0-(2.0d0/3.0d0)*ell*p20)
!            end if
!             dist=dsqrt(2-2*dsin(theta)*dsin(theta_surf)* &
!                  (dcos(phi)*dcos(phi_surf)+sin(phi)*sin(phi_surf)))
            dist= acos(cos(theta)*cos(theta_surf) + &
                       sin(theta)*sin(theta_surf)*cos(phi-phi_surf))*180.0d0/PI

            ! keep this point if it is closer to the receiver
            if(dist < distmin) then
              distmin = dist
              iglob_selected = iglob
              iface_selected = iface
              iselected = i
              jselected = j
              r_surf=dsqrt(dble(xstore(iglob_selected)**2)+ &
                     dble(ystore(iglob_selected)**2)+dble(zstore(iglob_selected)**2))
!              if(ELLIPTICITY) then
!                altitude_source(1) = r_surf-r_ellip
!              else
!                altitude_source(1) = r_surf-R_EARTH
!              end if
              altitude_rec(1) = r_surf-R_EARTH_SURF
            endif
          enddo
        enddo
      ! end of loop on all the elements on the free surface
      end do
      !  weighted mean at current point of topography elevation of the four closest nodes
      !  set distance to huge initial value
      distmin = HUGEVAL
      do j=jselected,jselected+1
        do i=iselected,iselected+1
          inode = 1
          do jadjust=0,1
            do iadjust= 0,1
              ispec = free_surface_ispec(iface_selected)
              igll = free_surface_ijk(1,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
              jgll = free_surface_ijk(2,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
              kgll = free_surface_ijk(3,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
              iglob = ibool(igll,jgll,kgll,ispec)
              r_surf=dsqrt(dble(xstore(iglob)**2)+dble(ystore(iglob)**2)+dble(zstore(iglob)**2))
              theta_surf=dacos(zstore(iglob)/r_surf)
!            sin(phi_surf)=ystore(iglob)/r_surf/dsin(theta_surf)
              phi_surf=dacos(dble(xstore(iglob))/(dble(r_surf)*dsin(theta_surf))*(1-1.e-7))
              if(ystore(iglob)<-1.e-14) then
                    phi_surf=-phi_surf+2*PI
              end if
!            if(ELLIPTICITY) then
!                dcost = dcos(theta_surf)
!                p20 = 0.5d0*(3.0d0*dcost*dcost-1.0d0)
!                call spline_evaluation(rspl,espl,espl2,nspl,R_EARTH,ell)
!                r_ellip = R_EARTH*(1.0d0-(2.0d0/3.0d0)*ell*p20)
!                elevation_node(inode) = r_surf-r_ellip
!              else
!                elevation_node(inode) = r_surf-R_EARTH
!              end if


              elevation_node(inode) = r_surf-R_EARTH_SURF
!              dist_node(inode) = dsqrt(2-2*dsin(theta)*dsin(theta_surf)* &
!                 (dcos(phi)*dcos(phi_surf)+sin(phi)*sin(phi_surf)))
              dist_node(inode)= acos(cos(theta)*cos(theta_surf) + &
                       sin(theta)*sin(theta_surf)*cos(phi-phi_surf))*180.0d0/PI
              inode = inode + 1
            end do
          end do
          dist = sum(dist_node)
          if(dist < distmin) then
            distmin = dist
            altitude_rec(1) = (dist_node(1)/dist)*elevation_node(1) + &
                       (dist_node(2)/dist)*elevation_node(2) + &
                       (dist_node(3)/dist)*elevation_node(3) + &
                       (dist_node(4)/dist)*elevation_node(4)
          endif
        end do
      end do
    end if
    !  MPI communications to determine the best slice
    distmin_ele(1)= distmin
    call gather_all_dp(distmin_ele,1,distmin_ele_all,1,NPROC)
    call gather_all_dp(altitude_rec,1,elevation_all,1,NPROC)
    if(myrank == 0) then
      iproc = minloc(distmin_ele_all)
      altitude_rec(1) = elevation_all(iproc(1))
    end if
    call bcast_all_dp(altitude_rec,1)
    elevation(irec) = altitude_rec(1)

    ! reset distance to huge initial value
    distmin=HUGEVAL

!     get the Cartesian components of n in the model: nu

    ! orientation consistent with the UTM projection
    ! X coordinate - East
    nu(1,1,irec) = 1.d0
    nu(1,2,irec) = 0.d0
    nu(1,3,irec) = 0.d0
    ! Y coordinate - North
    nu(2,1,irec) = 0.d0
    nu(2,2,irec) = 1.d0
    nu(2,3,irec) = 0.d0
    ! Z coordinate - Vertical
    nu(3,1,irec) = 0.d0
    nu(3,2,irec) = 0.d0
    nu(3,3,irec) = 1.d0

!    if(ELLIPTICITY) then
!      dcost = dcos(theta)
!      p20 = 0.5d0*(3.0d0*dcost*dcost-1.0d0)
!      call spline_evaluation(rspl,espl,espl2,nspl,R_EARTH,ell)
!      r_ellip  =  r_ellip*(1.0d0-(2.0d0/3.0d0)*ell*p20) 
!      r_target =  r_ellip_source + elevation(irec) - stbur(irec)
!    else
!      r_target =  R_EARTH + elevation(irec) - stbur(irec)
!    endif
    r_target(irec) = R_EARTH_SURF + elevation(irec) - stbur(irec)
    x_target(irec) = r_target(irec)*dsin(theta)*dcos(phi)
    y_target(irec) = r_target(irec)*dsin(theta)*dsin(phi)
    z_target(irec) = r_target(irec)*dcos(theta)

    !if (myrank == 0) write(IOVTK,*) x_target(irec), y_target(irec), z_target(irec)

    ! determines closest GLL point
    ispec_selected_rec(irec) = 0
    do ispec=1,NSPEC_AB

      ! define the interval in which we look for points
      if(FASTER_RECEIVERS_POINTS_ONLY) then
        imin = 1
        imax = NGLLX

        jmin = 1
        jmax = NGLLY

        kmin = 1
        kmax = NGLLZ

      else
        ! loop only on points inside the element
        ! exclude edges to ensure this point is not shared with other elements
        imin = 2
        imax = NGLLX - 1

        jmin = 2
        jmax = NGLLY - 1

        kmin = 2
        kmax = NGLLZ - 1
      endif

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            iglob = ibool(i,j,k,ispec)

            if (.not. RECEIVERS_CAN_BE_BURIED) then
              if ((.not. iglob_is_surface_external_mesh(iglob)) .or. (.not. ispec_is_surface_external_mesh(ispec))) then
                cycle
              endif
            endif

            dist = dsqrt((x_target(irec)-dble(xstore(iglob)))**2 &
                        +(y_target(irec)-dble(ystore(iglob)))**2 &
                        +(z_target(irec)-dble(zstore(iglob)))**2)

            ! keep this point if it is closer to the receiver
            if(dist < distmin) then
              distmin = dist
              ispec_selected_rec(irec) = ispec
!             if(DEBUG_COUPLING) print *,'ispec_selected_rec',ispec_selected_rec(irec),distmin,irec,myrank,x_target(irec),&
!                      y_target(irec),z_target(irec),elevation(irec),stbur(irec)
              ix_initial_guess(irec) = i
              iy_initial_guess(irec) = j
              iz_initial_guess(irec) = k

              xi_receiver(irec) = dble(ix_initial_guess(irec))
              eta_receiver(irec) = dble(iy_initial_guess(irec))
              gamma_receiver(irec) = dble(iz_initial_guess(irec))
              x_found(irec) = xstore(iglob)
              y_found(irec) = ystore(iglob)
              z_found(irec) = zstore(iglob)
            endif

          enddo
        enddo
      enddo
      ! compute final distance between asked and found (converted to km)
      final_distance(irec) = dsqrt((x_target(irec)-x_found(irec))**2 + &
        (y_target(irec)-y_found(irec))**2 + (z_target(irec)-z_found(irec))**2)

    ! end of loop on all the spectral elements in current slice
    enddo
    if (ispec_selected_rec(irec) == 0) then
      final_distance(irec) = HUGEVAL
    endif
    ! record three components for each station
   do iorientation = 1,3

!     North
       if(iorientation == 1) then
        stazi = 0.d0
        stdip = 0.d0
!     East
       else if(iorientation == 2) then
         stazi = 90.d0
         stdip = 0.d0
!     Vertical
       else if(iorientation == 3) then
         stazi = 0.d0
         stdip = - 90.d0
       else
         call exit_MPI(myrank,'incorrect orientation')
       endif
!     get the orientation of the seismometer
      thetan=(90.0d0+stdip)*PI/180.0d0
      phin=stazi*PI/180.0d0

! we use the same convention as in Harvard normal modes for the orientation
      !     vertical component
      n(1) = cos(thetan)
      !     N-S component
      n(2) = - sin(thetan)*cos(phin)
      !     E-W component
      n(3) = sin(thetan)*sin(phin)

!     get the Cartesian components of n in the model: nu
      sint = sin(theta)
      cost = cos(theta)
      sinp = sin(phi)
      cosp = cos(phi)
      

      ! build rotation matrice nu for seismograms
      nu(iorientation,1,irec) = n(1)*sint*cosp+n(2)*cost*cosp-n(3)*sinp
      nu(iorientation,2,irec) = n(1)*sint*sinp+n(2)*cost*sinp+n(3)*cosp
      nu(iorientation,3,irec) = n(1)*cost-n(2)*sint
      if(DEBUG_COUPLING) print *,'ENZ',irec,iorientation,nu(iorientation,:,irec)

   !end of loop on iorientation
   enddo


  ! end of loop on all the stations
  enddo

  ! close receiver file
  close(1)

! ****************************************
! find the best (xi,eta,gamma) for each receiver
! ****************************************

  if(.not. FASTER_RECEIVERS_POINTS_ONLY) then

    ! loop on all the receivers to iterate in that slice
    do irec = 1,nrec
      ispec_iterate = ispec_selected_rec(irec)

      ! use initial guess in xi and eta
      xi = xigll(ix_initial_guess(irec))
      eta = yigll(iy_initial_guess(irec))
      gamma = zigll(iz_initial_guess(irec))

      ! define coordinates of the control points of the element
      do ia=1,NGNOD
        iax = 0
        iay = 0
        iaz = 0
        if(iaddx(ia) == 0) then
          iax = 1
        else if(iaddx(ia) == 1) then
          iax = (NGLLX+1)/2
        else if(iaddx(ia) == 2) then
          iax = NGLLX
        else
          call exit_MPI(myrank,'incorrect value of iaddx')
        endif

        if(iaddy(ia) == 0) then
          iay = 1
        else if(iaddy(ia) == 1) then
          iay = (NGLLY+1)/2
        else if(iaddy(ia) == 2) then
          iay = NGLLY
        else
          call exit_MPI(myrank,'incorrect value of iaddy')
        endif

        if(iaddz(ia) == 0) then
          iaz = 1
        else if(iaddz(ia) == 1) then
          iaz = (NGLLZ+1)/2
        else if(iaddz(ia) == 2) then
          iaz = NGLLZ
        else
          call exit_MPI(myrank,'incorrect value of iaddz')
        endif

        iglob = ibool(iax,iay,iaz,ispec_iterate)
        xelm(ia) = dble(xstore(iglob))
        yelm(ia) = dble(ystore(iglob))
        zelm(ia) = dble(zstore(iglob))

      enddo
     
      ! iterate to solve the non linear system
      do iter_loop = 1,NUM_ITER

        ! impose receiver exactly at the surface
        !    gamma = 1.d0

        ! recompute jacobian for the new point
        call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,NGNOD)

        ! compute distance to target location
        dx = - (x - x_target(irec))
        dy = - (y - y_target(irec))
        dz = - (z - z_target(irec))

        ! compute increments
        ! gamma does not change since we know the receiver is exactly on the surface
        dxi  = xix*dx + xiy*dy + xiz*dz
        deta = etax*dx + etay*dy + etaz*dz
        dgamma = gammax*dx + gammay*dy + gammaz*dz
        ! update values
        xi = xi + dxi
        eta = eta + deta
        gamma = gamma + dgamma

        ! impose that we stay in that element
        ! (useful if user gives a receiver outside the mesh for instance)
        ! we can go slightly outside the [1,1] segment since with finite elements
        ! the polynomial solution is defined everywhere
        ! this can be useful for convergence of itertive scheme with distorted elements
        if (xi > 1.10d0) xi = 1.10d0
        if (xi < -1.10d0) xi = -1.10d0
        if (eta > 1.10d0) eta = 1.10d0
        if (eta < -1.10d0) eta = -1.10d0
        if (gamma > 1.10d0) gamma = 1.10d0
        if (gamma < -1.10d0) gamma = -1.10d0

      ! end of non linear iterations
      enddo

      ! impose receiver exactly at the surface after final iteration
      !  gamma = 1.d0

      ! compute final coordinates of point found
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,NGNOD)

      ! store xi,eta and x,y,z of point found
      xi_receiver(irec) = xi
      eta_receiver(irec) = eta
      gamma_receiver(irec) = gamma
      x_found(irec) = x
      y_found(irec) = y
      z_found(irec) = z

      ! compute final distance between asked and found (converted to km)
      final_distance(irec) = dsqrt((x_target(irec)-x_found(irec))**2 + &
        (y_target(irec)-y_found(irec))**2 + (z_target(irec)-z_found(irec))**2)

    enddo

  endif ! of if (.not. FASTER_RECEIVERS_POINTS_ONLY)

  ! synchronize all the processes to make sure all the estimates are available
  call synchronize_all()

  ! for MPI version, gather information from all the nodes
  ispec_selected_rec_all(:,:) = -1
  call gather_all_i(ispec_selected_rec,nrec,ispec_selected_rec_all,nrec,NPROC)
  call gather_all_dp(xi_receiver,nrec,xi_receiver_all,nrec,NPROC)
  call gather_all_dp(eta_receiver,nrec,eta_receiver_all,nrec,NPROC)
  call gather_all_dp(gamma_receiver,nrec,gamma_receiver_all,nrec,NPROC)
  call gather_all_dp(final_distance,nrec,final_distance_all,nrec,NPROC)
  call gather_all_dp(x_found,nrec,x_found_all,nrec,NPROC)
  call gather_all_dp(y_found,nrec,y_found_all,nrec,NPROC)
  call gather_all_dp(z_found,nrec,z_found_all,nrec,NPROC)
  call gather_all_dp(nu,3*3*nrec,nu_all,3*3*nrec,NPROC)

  ! this is executed by main process only
  if(myrank == 0) then

    ! check that the gather operation went well
    if(any(ispec_selected_rec_all(:,:) == -1)) call exit_MPI(myrank,'gather operation failed for receivers')

    ! MPI loop on all the results to determine the best slice
    islice_selected_rec(:) = -1
    do irec = 1,nrec
    distmin = HUGEVAL
    do iprocloop = 0,NPROC-1
      if(final_distance_all(irec,iprocloop) < distmin) then
        distmin = final_distance_all(irec,iprocloop)
        islice_selected_rec(irec) = iprocloop
        ispec_selected_rec(irec) = ispec_selected_rec_all(irec,iprocloop)
        xi_receiver(irec) = xi_receiver_all(irec,iprocloop)
        eta_receiver(irec) = eta_receiver_all(irec,iprocloop)
        gamma_receiver(irec) = gamma_receiver_all(irec,iprocloop)
        x_found(irec) = x_found_all(irec,iprocloop)
        y_found(irec) = y_found_all(irec,iprocloop)
        z_found(irec) = z_found_all(irec,iprocloop)
        nu(:,:,irec) = nu_all(:,:,irec,iprocloop)
      endif
    enddo
    final_distance(irec) = distmin
    enddo
    nrec_found = 0
    do irec=1,nrec
    
      if(final_distance(irec) == HUGEVAL) call exit_MPI(myrank,'error locating receiver')

      if(final_distance(irec) > 3000.d0) then
        write(IMAIN,*) 'station # ',irec,'     ',station_name(irec),network_name(irec)
        write(IMAIN,*) '*****************************************************************'
        write(IMAIN,*) '***** WARNING: receiver is located outside the mesh, therefore excluded *****'
      else
        nrec_found = nrec_found + 1

        islice_selected_rec_found(nrec_found) = islice_selected_rec(irec)
        ispec_selected_rec_found(nrec_found) = ispec_selected_rec(irec)
        xi_receiver_found(nrec_found) = xi_receiver(irec)
        eta_receiver_found(nrec_found) = eta_receiver(irec)
        gamma_receiver_found(nrec_found) = gamma_receiver(irec)
        station_name_found(nrec_found) = station_name(irec)
        network_name_found(nrec_found) = network_name(irec)
        stlat_found(nrec_found) = stlat(irec)
        stlon_found(nrec_found) = stlon(irec)
        stele_found(nrec_found) = stele(irec)
        stbur_found(nrec_found) = stbur(irec)
        nu_found(:,:,nrec_found) = nu(:,:,irec)
        epidist_found(nrec_found) = epidist(irec)
      endif
    end do
    ! compute maximal distance for all the receivers
    final_distance_max = maxval(final_distance(:))
    ! display maximum error for all the receivers
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of all the receivers: ',sngl(final_distance_max),' km'
    if(nrec.le.0) call exit_MPI(myrank,'At least one station should be in simulated region')
    nrec = nrec_found
    islice_selected_rec(1:nrec) = islice_selected_rec_found(1:nrec)
    ispec_selected_rec(1:nrec) = ispec_selected_rec_found(1:nrec)
    xi_receiver(1:nrec) = xi_receiver_found(1:nrec)
    eta_receiver(1:nrec) = eta_receiver_found(1:nrec)
    gamma_receiver(1:nrec) = gamma_receiver_found(1:nrec)
    station_name(1:nrec) = station_name_found(1:nrec)
    network_name(1:nrec) = network_name_found(1:nrec)
    stlat(1:nrec) = stlat_found(1:nrec)
    stlon(1:nrec) = stlon_found(1:nrec)
    stele(1:nrec) = stele_found(1:nrec)
    stbur(1:nrec) = stbur_found(1:nrec)
    nu(:,:,1:nrec) = nu_found(:,:,1:nrec)
    epidist(1:nrec) = epidist_found(1:nrec)


    do irec=1,nrec

      write(IMAIN,*)
      write(IMAIN,*) 'station # ',irec,'    ',station_name(irec),network_name(irec)

      if(final_distance(irec) == HUGEVAL) call exit_MPI(myrank,'error locating receiver')

      write(IMAIN,*) '     original latitude: ',sngl(stlat(irec))
      write(IMAIN,*) '     original longitude: ',sngl(stlon(irec))
      write(IMAIN,*) '     original depth: ',sngl(stbur(irec)),' m'

      write(IMAIN,*) '     target x, y, z, r: ',sngl(x_target(irec)),sngl(y_target(irec)), &
                           sngl(z_target(irec)),sngl(r_target(irec))
      write(IMAIN,*) '     distance: ',epidist(irec),' deg'

      write(IMAIN,*) '     closest estimate found: ',sngl(final_distance(irec)),' m away'
      write(IMAIN,*) '     in slice ',islice_selected_rec(irec),' in element ',ispec_selected_rec(irec)
      if(FASTER_RECEIVERS_POINTS_ONLY) then
        write(IMAIN,*) '     in point i,j,k = ',nint(xi_receiver(irec)),nint(eta_receiver(irec)),nint(gamma_receiver(irec))
        write(IMAIN,*) '     nu1 = ',nu(1,:,irec)
        write(IMAIN,*) '     nu2 = ',nu(2,:,irec)
        write(IMAIN,*) '     nu3 = ',nu(3,:,irec)
      else
        write(IMAIN,*) '     at coordinates: '
        write(IMAIN,*) '       xi    = ',xi_receiver(irec)
        write(IMAIN,*) '       eta   = ',eta_receiver(irec)
        write(IMAIN,*) '       gamma = ',gamma_receiver(irec)
      endif
      write(IMAIN,*) '         x: ',x_found(irec)
      write(IMAIN,*) '         y: ',y_found(irec)
      write(IMAIN,*) '         z: ',z_found(irec)
      write(IMAIN,*)


      ! add warning if estimate is poor
      ! (usually means receiver outside the mesh given by the user)
      if(final_distance(irec) > 3000.d0) then
        write(IMAIN,*) '*******************************************************'
        write(IMAIN,*) '***** WARNING: receiver location estimate is poor *****'
        write(IMAIN,*) '*******************************************************'
      endif

      write(IMAIN,*)

    enddo

    ! compute maximal distance for all the receivers
    final_distance_max = maxval(final_distance(1:nrec))

    ! display maximum error for all the receivers
    write(IMAIN,*) 'maximum error in location of all the receivers: ',sngl(final_distance_max),' m'

    ! add warning if estimate is poor
    ! (usually means receiver outside the mesh given by the user)
    if(final_distance_max > 1000.d0) then
      write(IMAIN,*)
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '***** WARNING: at least one receiver is poorly located *****'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
    endif

    ! get the base pathname for output files
    !call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)))

    ! write the list of stations and associated epicentral distance
    open(unit=27,file=trim(OUTPUT_FILES)//'/output_list_stations.txt',status='unknown')
    do irec=1,nrec
      write(27,*) station_name(irec),'.',network_name(irec),' : ',epidist(irec),' deg horizontal distance'
    enddo
    close(27)

    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of receiver detection - done'
    write(IMAIN,*)

  endif    ! end of section executed by main process only

  ! main process broadcasts the results to all the slices
  call bcast_all_i(islice_selected_rec,nrec)
  call bcast_all_i(ispec_selected_rec,nrec)
  call bcast_all_dp(xi_receiver,nrec)
  call bcast_all_dp(eta_receiver,nrec)
  call bcast_all_dp(gamma_receiver,nrec)
  ! synchronize all the processes to make sure everybody has finished
  call synchronize_all()

  ! deallocate arrays
  deallocate(stlat)
  deallocate(stlon)
  deallocate(stele)
  deallocate(stbur)
  deallocate(epidist)
  deallocate(ix_initial_guess)
  deallocate(iy_initial_guess)
  deallocate(iz_initial_guess)
  deallocate(r_target)
  deallocate(x_target)
  deallocate(y_target)
  deallocate(z_target)
  deallocate(x_found)
  deallocate(y_found)
  deallocate(z_found)
  deallocate(final_distance)
  deallocate(islice_selected_rec_found)
  deallocate(ispec_selected_rec_found)
  deallocate(xi_receiver_found)
  deallocate(eta_receiver_found)
  deallocate(gamma_receiver_found)
  deallocate(nu_found)
  deallocate(station_name_found)
  deallocate(network_name_found)
  deallocate(stlat_found)
  deallocate(stlon_found)
  deallocate(stele_found)
  deallocate(stbur_found)
  deallocate(epidist_found)
  deallocate(ispec_selected_rec_all)
  deallocate(xi_receiver_all)
  deallocate(eta_receiver_all)
  deallocate(gamma_receiver_all)
  deallocate(x_found_all)
  deallocate(y_found_all)
  deallocate(z_found_all)
  deallocate(final_distance_all)

  end subroutine locate_receivers_cubed

!=====================================================================

