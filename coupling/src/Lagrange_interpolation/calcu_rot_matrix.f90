subroutine calcu_rotmatrix_idGreen ()
   use convolution_par
   use constants
   use coupling_SEM_DSM_par,only:nproc,myrank,npack,fft_disp_bound_new,fft_traction_bound_new,&
                          normal_vector,coord_bound,c,jacobian,NPOINTS,id_depth,is_elastic,&
                          is_acoustic
   implicit none

!other variables
   integer ::ista,ipoint,icomp
   double precision,dimension(ncomp) ::vector_z_bound,vector_r_bound,vector_t_bound
   double precision,dimension(ncomp) ::vector_z_sta,vector_r_sta,vector_t_sta
   double precision ::theta,phi,radius,depth
   double precision ::cos_gcarc,gcarc
   integer ::idep,igcarc
   double precision ::modu_temp


   do ista=1,nstation
    theta=pi/2.d0-station_lat(ista)*degtorad
    phi=station_lon(ista)*degtorad
    station_coord(1,ista)=dsin(theta)*dcos(phi)
    station_coord(2,ista)=dsin(theta)*dsin(phi)
    station_coord(3,ista)=dcos(theta)
    vector_z_sta(:)=station_coord(:,ista)
    do ipoint=1,npoints
      radius=dsqrt(dot_product(coord_bound(:,ipoint),coord_bound(:,ipoint)))
      vector_z_bound(:)=coord_bound(:,ipoint)
      depth=r_earth-radius
      idep=int((depth-min_dep)/ddep+1)

!To be removed
!      if(idep.gt.-8.and.idep.lt.1) idep=1

      do icomp=1,ncomp
         vector_z_bound(icomp)=vector_z_bound(icomp)/radius
      end do
!      statopoint(:)=station_coord(:,ista)-vector_r
      vector_t_bound(1)= station_coord(2,ista)*vector_z_bound(3) - &
                          station_coord(3,ista)*vector_z_bound(2)
      vector_t_bound(2)= -station_coord(1,ista)*vector_z_bound(3)+&
                          station_coord(3,ista)*vector_z_bound(1)
      vector_t_bound(3)=  station_coord(1,ista)*vector_z_bound(2)- &
                          station_coord(2,ista)*vector_z_bound(1)
      modu_temp=dsqrt(dot_product(vector_t_bound,vector_t_bound))
      do icomp=1,ncomp
        vector_t_bound(icomp)=vector_t_bound(icomp)/modu_temp
      end do
      vector_t_sta(:)=vector_t_bound(:)

      vector_r_bound(1)=vector_t_bound(2)*vector_z_bound(3)- &
                        vector_t_bound(3)*vector_z_bound(2)
      vector_r_bound(2)=-vector_t_bound(1)*vector_z_bound(3)+ &
                        vector_t_bound(3)*vector_z_bound(1)
      vector_r_bound(3)=vector_t_bound(1)*vector_z_bound(2)- &
                        vector_t_bound(2)*vector_z_bound(1)
      vector_r_sta(1)=vector_t_sta(2)*vector_z_sta(3)- &
                        vector_t_sta(3)*vector_z_sta(2)
      vector_r_sta(2)=-vector_t_sta(1)*vector_z_sta(3)+ &
                        vector_t_sta(3)*vector_z_sta(1)
      vector_r_sta(3)=vector_t_sta(1)*vector_z_sta(2)- &
                        vector_t_sta(2)*vector_z_sta(1)

      rot_matrix_bound(1,:,ipoint,ista)=vector_z_bound(:)
      rot_matrix_bound(2,:,ipoint,ista)=vector_r_bound(:)
      rot_matrix_bound(3,:,ipoint,ista)=vector_t_bound(:)

      rot_matrix_station(1,:,ipoint,ista)=vector_z_sta(:)
      rot_matrix_station(2,:,ipoint,ista)=vector_r_sta(:)
      rot_matrix_station(3,:,ipoint,ista)=vector_t_sta(:)


      cos_gcarc=dot_product(vector_z_sta,vector_z_bound)
      gcarc=dacos(cos_gcarc)*radtodeg
      igcarc=int((dabs(gcarc)-min_theta)/dtheta+1)


!to be fixed
      if(is_acoustic(ipoint).and. &
         (id_depth(ipoint).lt.1.or.id_depth(ipoint).gt.ndep_acous)) then
      !   call exit_MPI('Error, idep should be between 1 and ndep_acous')
      end if

      if(is_elastic(ipoint).and.&
         (id_depth(ipoint).lt.1.or.id_depth(ipoint).gt.ndep_elas)) then
         call exit_MPI('Error, idep should be between 1 and ndep_elas')
      end if

      if(igcarc.le.nfit_Green/2.or.igcarc.ge.ntheta-nfit_Green/2) then
         call exit_MPI('Error, igcarc should be between 1 and ntheta-1')
      end if

      id_Green(ipoint,ista)=ntheta*(id_depth(ipoint)-1)+igcarc
      distance(ipoint,ista)=gcarc
      idistance(ipoint,ista)=igcarc
!      if(myrank.eq.0) print *,'id_Green',gcarc,igcarc,igcarc_remained(ipoint,ista),id_Green(ipoint,ista)

    end do
   end do


end subroutine calcu_rotmatrix_idGreen
