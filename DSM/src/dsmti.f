	program dsmti
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computation of synthetic seismograms for spherical symmetric
c TI media by using the Direct Solution Method.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c constant
        implicit none
	integer maxngrid_r,maxlmax,nl_check_amp
	integer max_nstation,max_ndep,maxn_structure_zone,maxnfreq
        integer max_ntheta
	parameter ( maxngrid_r = 100000 )
	parameter ( maxlmax = 33000 )
c WENBO
        parameter (max_ndep = 400 )
        parameter (max_ntheta=2400)
	parameter ( max_nstation = max_ndep*max_ntheta )
c WENBO --- nl_check_amp need to be tested. 
c       --- Initial tests shows that it would be nice to set it larger than 400
        parameter ( nl_check_amp=500)
	parameter ( maxn_structure_zone = 170 )
	parameter ( maxnfreq = 200 )

	real*8 pi
	parameter ( pi=3.1415926535897932d0 )
c variables for time series
	integer n_frequency,i_frequency
	real*8 time_series_length,omega_imag
c variables for numerical grids
	integer ngrid_r0,lmin0,lmax0
	integer ngrid_r_estimated,ngrid_r,lmin,lmax,l
        integer lmax_computed,lmax_allfreq,ierr
	integer n_frequency_band,i_frequency_band
	integer i_frequency_min,i_frequency_max
	real*8  grid_r(maxngrid_r)
        real*8  r_CMB,r_ICB
        integer ir_CMB,ir_ICB
c WENBO
        integer idim_ir_sph(maxngrid_r),idim_ir_tor(maxngrid_r)
        integer idim_ir_sph0(maxngrid_r)
	real*8 grid_rho(2,maxngrid_r)
	real*8 grid_Ak(2,maxngrid_r),grid_Am(2,maxngrid_r)
	real*8 grid_Ck(2,maxngrid_r),grid_Cm(2,maxngrid_r)
	real*8 grid_L(2,maxngrid_r),grid_N(2,maxngrid_r)
	real*8 grid_Fk(2,maxngrid_r),grid_Fm(2,maxngrid_r)
	real*8 grid_kappa(2,maxngrid_r)
	real*8 grid_mu(2,maxngrid_r)
	real*8 grid_qkappa(maxngrid_r)
	real*8 grid_qmu(maxngrid_r)
c variables for the structure
	integer n_structure_zone
	integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
	real*8 rmin_structure_zone(maxn_structure_zone)
	real*8 rmax_structure_zone(maxn_structure_zone)
	real*8 rho_structure_zone(4,maxn_structure_zone)
	real*8 vpv_structure_zone(4,maxn_structure_zone)
	real*8 vph_structure_zone(4,maxn_structure_zone)
	real*8 vsv_structure_zone(4,maxn_structure_zone)
	real*8 vsh_structure_zone(4,maxn_structure_zone)
	real*8 eta_structure_zone(4,maxn_structure_zone)
	real*8 qkappa_structure_zone(maxn_structure_zone)
	real*8 qmu_structure_zone(maxn_structure_zone)
!WENBO  
        integer fluid_thiszone(maxn_structure_zone)
c variables for a source
	integer igrid_rs
	integer idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0
	real*8 source_r,source_mt(3,3)
c WENBO 
        integer source_type
        real*8 fr,ftheta,fphi
	real*8 source_depth,source_lat,source_lon
	real*8 grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
	real*8 grid_qkappas,grid_qmus
c variables for stations
	integer n_station_solid,n_station_fluid
c WENBO
c        real*8  min_dep_solid,min_dep_fluid,ddep
        real*8  depth_solid(max_ndep),depth_fluid(max_ndep)
        integer izone_idep_solid(max_ndep)
        integer izone_idep_fluid(max_ndep)
        real*8  depth_for_stress(3,max_ndep)
        real*8  depth_for_pressure(3,max_ndep)
        real*8  r_top_stress(3),r_bot_stress(3)
        real*8  r_top_pressure(3),r_bot_pressure(3)
        integer ir_top_stress(3),ir_bot_stress(3)
        integer ir_top_pressure(3),ir_bot_pressure(3)
        real*8  depth_tolerence
        integer ndep_solid,ndep_fluid,idep_forcut
        integer ndist_solid,ndist_fluid
        real*8  dist_solid(max_ntheta),dist_fluid(max_ntheta)
        integer top_fluid,bot_fluid
c        real*8  min_theta,dtheta
c        integer ntheta
        integer save_velo
        integer ir_dep_solid(max_ndep),ir_dep_fluid(max_ndep)
        integer ir_for_stress(3,max_ndep)
        integer ir_for_pressure(3,max_ndep)
        real*8  r_sta_solid(max_ndep),r_sta_fluid(max_ndep)
c	integer idim_station_sph,idim_station_tor,
c    &	        idim_station_sph0,idim_station_tor0
        integer idim_sph_forcut,minidim_sph_sta,
     &          idim_sph_forcut0,minidim_sph_sta0,
     &          idim_tor_forcut,minidim_tor_sta,
     &          idim_tor_forcut0,minidim_tor_sta0
        integer min_ir_dep
c	real*8 station_lat(max_nstation),station_lon(max_nstation)
c	real*8 station_theta(max_nstation),station_phi(max_nstation)
c	complex*16 vecsph_sph1(3,0:maxlmax,-2:2,max_nstation)
c	complex*16 vecsph_sph2(3,0:maxlmax,-2:2,max_nstation)
c	complex*16 vecsph_tor(3,0:maxlmax,-2:2,max_nstation)
	complex*16 station_disp_solid(3,max_nstation),
     &             station_disp_fluid(3,max_nstation)
        complex*16 coef_c(0:maxlmax,-2:2,2,3)
        complex*16 coef_dcdr(0:maxlmax,-2:2,2,3)
        real*8 Atop,Ctop,Ftop,Ltop,Ntop,
     &         Abot,Cbot,Fbot,Lbot,Nbot
c WENBO
        complex*16 epsilon11(max_nstation),
     &            epsilon22(max_nstation),
     &            epsilon33(max_nstation),
     &            epsilon12(max_nstation),
     &            epsilon13(max_nstation),
     &            epsilon23(max_nstation)
        complex*16 pressure(max_nstation)
c	character*80 sac_file(max_nstation)
c WENBO

c variables for matrix elements
	real*8 submatrix_I0(4,maxngrid_r)
	real*8 submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
	real*8 submatrix_I2(4,maxngrid_r)
	real*8 submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
	real*8 submatrix_I4(4,maxngrid_r)
	real*8 submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
	real*8 submatrix_I6(4,maxngrid_r)
	real*8 submatrix_I7(4,maxngrid_r)
	real*8 submatrix_I3k_mod(6,maxngrid_r)
	real*8 submatrix_I3m_mod(6,maxngrid_r)
	real*8 submatrix_I4_mod(6,maxngrid_r)
c variables for wavefield
	integer i_significance,idim0,init_npos_sph,init_npos_tor
	real*8 amp_max
	complex*16 omega
	complex*16 whole_matrix_sph(4,2*maxngrid_r)
	complex*16 whole_matrix_tor(2,maxngrid_r)
	complex*16 whole_matrix_dr_sph(2*maxngrid_r)
	complex*16 whole_matrix_dr_tor(maxngrid_r)
	complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
	complex*16 whole_vector_tor(maxngrid_r,-2:2)
	complex*16 work_vector(2*maxngrid_r)
	real work_time(32*maxnfreq)
c	complex*16 work_spc(16*maxnfreq)
c	equivalence ( work_vector,work_spc )
c MPI variables
	include 'mpicom.h'
	integer info,irank,nsnd
c functions
	omega(time_series_length,i_frequency,omega_imag)
     &	  = dcmplx( 2.d0*pi*dble(i_frequency)/dble(time_series_length),
     &	            -omega_imag )
c	ngrid_r_estimated
c     &	  (i_frequency_band,n_frequency_band,ngrid_r0)
c     &	  = int(
c     &	      dnint (
c     &	        ngrid_r0
c     &	          * dble(i_frequency_band) / dble(n_frequency_band)
c     &	      )
c     &	    )
        ngrid_r_estimated
     &    (i_frequency_band,n_frequency_band,ngrid_r0)
     &    = int(
     &        dnint (
     &          ngrid_r0
     &            * dble(n_frequency_band) / dble(n_frequency_band)
     &        )
     &      )

	lmin
     &	  (i_frequency_band,n_frequency_band,lmin0 )
     &    = lmin0
	lmax
     &	  (i_frequency_band,n_frequency_band,lmax0 )
     &	  = int(
     &	      dnint (
     &	        lmax0
     &	          * dble(i_frequency_band) / dble(n_frequency_band)
     &	      )
     &	    )
	i_frequency_min
     &	  (i_frequency_band,n_frequency_band,n_frequency)
     &	  = int(
     &	      dnint (
     &	        n_frequency
     &	          * dble(i_frequency_band-1) / dble(n_frequency_band)
     &	      )
     &	    ) + 1 - int(1.d0/dble(i_frequency_band))
	i_frequency_max
     &	  (i_frequency_band,n_frequency_band,n_frequency)
     &	  = int(
     &	      dnint (
     &	        n_frequency
     &	          * dble(i_frequency_band) / dble(n_frequency_band)
     &	      )
     &	    )
c
	data n_frequency_band / 4 /
c **********************************************************************
c inputing parameters
c **********************************************************************
c initialize MPI
	call minit

c
c WENBO
	if ( rank.eq.0 ) then
	  call param_input
     &	  ( maxngrid_r,maxlmax,max_nstation,maxn_structure_zone,
     &	    time_series_length,n_frequency,omega_imag,
     &	    max_ndep,max_ntheta,ngrid_r0,lmin0,lmax0,
     &	    n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &	    rho_structure_zone,vpv_structure_zone,vph_structure_zone,
     &	    vsv_structure_zone,vsh_structure_zone,eta_structure_zone,
     &	    qkappa_structure_zone,qmu_structure_zone,
     &      fluid_thiszone,r_CMB,r_ICB,
     &	    source_r,source_type,source_mt,fr,ftheta,fphi,
     &      source_depth,source_lat,source_lon,
     &      depth_solid,depth_fluid,ndist_solid,
     &      ndist_fluid,dist_solid,dist_fluid,
     &      ndep_solid,ndep_fluid,depth_for_stress,
     &      depth_for_pressure,izone_idep_solid,izone_idep_fluid,
     &      depth_tolerence,top_fluid,bot_fluid,
     &      idep_forcut,save_velo,n_station_solid,
     &      n_station_fluid)
	  call msndi
	  call mputi( n_frequency,1,info )
	  call mputi( ngrid_r0,1,info )
	  call mputi( lmin0,1,info )
	  call mputi( lmax0,1,info )
	  call mputi( n_structure_zone,1,info )
	  call mputi( n_station_solid,1,info )
          call mputi( n_station_fluid,1,info )
          call mputi( ndep_solid,1,info )
          call mputi( ndep_fluid,1,info )
          call mputi( izone_idep_solid,max_ndep,info)
          call mputi( izone_idep_fluid,max_ndep,info)
          call mputi( idep_forcut,1,info )
          call mputi( ndist_solid,1,info )
          call mputi( ndist_fluid,1,info )
c WENBO
          call mputi( source_type,1,info)
          call mputi( save_velo,1,info)
          call mputi( fluid_thiszone,maxn_structure_zone,info )
          call mputi( top_fluid,1,info)
          call mputi( bot_fluid,1,info)

c WENBO
	  call mputd( time_series_length,1,info )
	  call mputd( omega_imag,1,info )
	  call mputd( rmin_structure_zone,maxn_structure_zone,info )
	  call mputd( rmax_structure_zone,maxn_structure_zone,info )
	  call mputd( rho_structure_zone,4*maxn_structure_zone,info )
	  call mputd( vpv_structure_zone,4*maxn_structure_zone,info )
	  call mputd( vph_structure_zone,4*maxn_structure_zone,info )
	  call mputd( vsv_structure_zone,4*maxn_structure_zone,info )
	  call mputd( vsh_structure_zone,4*maxn_structure_zone,info )
	  call mputd( eta_structure_zone,4*maxn_structure_zone,info )
	  call mputd( qkappa_structure_zone,maxn_structure_zone,info )
	  call mputd( qmu_structure_zone,maxn_structure_zone,info )
          call mputd( r_CMB,1,info )
          call mputd( r_ICB,1,info )
	  call mputd( source_r,1,info )
	  call mputd( source_mt,9,info )
c WENBO
          call mputd( fr,1,info)
          call mputd( ftheta,1,info)
          call mputd( fphi,1,info)
          call mputd( depth_tolerence,1,info)
          call mputd( depth_solid,max_ndep,info)
          call mputd( depth_fluid,max_ndep,info)
          call mputd( dist_solid,max_ntheta,info)
          call mputd( dist_fluid,max_ntheta,info)
          call mputd( depth_for_stress,3*max_ndep,info)
          call mputd( depth_for_pressure,3*max_ndep,info)

c WENBO
c	  call mputd( station_theta,max_nstation,info )
c	  call mputd( station_phi,max_nstation,info )
	  call msend( -1,10,info )
	else
	  call mrecv( 10,0,info )
	  call mgeti( n_frequency,1,info )
	  call mgeti( ngrid_r0,1,info )
	  call mgeti( lmin0,1,info )
	  call mgeti( lmax0,1,info )
	  call mgeti( n_structure_zone,1,info )
	  call mgeti( n_station_solid,1,info )
          call mgeti( n_station_fluid,1,info )
          call mgeti( ndep_solid,1,info)
          call mgeti( ndep_fluid,1,info)
          call mgeti( izone_idep_solid,max_ndep,info)
          call mgeti( izone_idep_fluid,max_ndep,info)
          call mgeti( idep_forcut,1,info )
          call mgeti( ndist_solid,1,info )
          call mgeti( ndist_fluid,1,info )

c WENBO
          call mgeti( source_type,1,info)
          call mgeti( save_velo,1,info)
          call mgeti( fluid_thiszone,maxn_structure_zone,info )

          call mgeti( top_fluid,1,info)
          call mgeti( bot_fluid,1,info)
c WENBO
	  call mgetd( time_series_length,1,info )
	  call mgetd( omega_imag,1,info )
	  call mgetd( rmin_structure_zone,maxn_structure_zone,info )
	  call mgetd( rmax_structure_zone,maxn_structure_zone,info )
	  call mgetd( rho_structure_zone,4*maxn_structure_zone,info )
	  call mgetd( vpv_structure_zone,4*maxn_structure_zone,info )
	  call mgetd( vph_structure_zone,4*maxn_structure_zone,info )
	  call mgetd( vsv_structure_zone,4*maxn_structure_zone,info )
	  call mgetd( vsh_structure_zone,4*maxn_structure_zone,info )
	  call mgetd( eta_structure_zone,4*maxn_structure_zone,info )
	  call mgetd( qkappa_structure_zone,maxn_structure_zone,info )
	  call mgetd( qmu_structure_zone,maxn_structure_zone,info )
          call mgetd( r_CMB,1,info)
          call mgetd( r_ICB,1,info)
	  call mgetd( source_r,1,info )
	  call mgetd( source_mt,9,info )
c WENBO
          call mgetd( fr,1,info)
          call mgetd( ftheta,1,info)
          call mgetd( fphi,1,info)
c WENBO
          call mgetd( depth_tolerence,1,info)
          call mgetd( depth_solid,max_ndep,info)
          call mgetd( depth_fluid,max_ndep,info)
          call mgetd( dist_solid,max_ntheta,info)
          call mgetd( dist_fluid,max_ntheta,info)
          call mgetd( depth_for_stress,3*max_ndep,info)
          call mgetd( depth_for_pressure,3*max_ndep,info)

c	  call mgetd( station_theta,max_nstation,info )
c	  call mgetd( station_phi,max_nstation,info )
	endif
c
c	call comp_vecsph_station
c     &	  ( maxlmax,lmax0,n_station,station_theta,station_phi,
c     &	    vecsph_sph1,vecsph_sph2,vecsph_tor )

c        if(rank.eq.0) call save_par(n_frequency,rank,   
c     &                 n_station_solid,n_station_fluid,
c     &                 ndist_solid,ndist_fluid,
c     &                 ndep_solid,ndep_fluid )

c
	do 1000 i_frequency_band=1,n_frequency_band
c **********************************************************************
c generating numerical grids
c **********************************************************************
           call grid_generation
     &      (maxngrid_r,ngrid_r_estimated
     &          (i_frequency_band,n_frequency_band,ngrid_r0),
     &      n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &      vsv_structure_zone,vsh_structure_zone,
     &      fluid_thiszone,r_CMB,r_ICB,ir_CMB,ir_ICB,
     &      grid_r,source_r,ngrid_r,depth_solid,depth_fluid,
     &      ndep_solid,ndep_fluid,max_ndep,depth_for_stress,
     &      depth_for_pressure,ir_for_stress,ir_for_pressure,
     &      ir_dep_solid,ir_dep_fluid,r_sta_solid,r_sta_fluid,
     &      depth_tolerence,izone_idep_solid,izone_idep_fluid,
     &      r_top_stress,r_bot_stress,
     &      r_top_pressure,r_bot_pressure,
     &      ir_top_stress,ir_bot_stress,ir_top_pressure,
     &      ir_bot_pressure,top_fluid,bot_fluid)

          print *,'gene_grid done, rank=',rank,',ifreq_band=',
     &        i_frequency_band
	  call assign_structure
     &	    ( ngrid_r,grid_r,
     &	      n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &	      rho_structure_zone,
     &	      vpv_structure_zone,vph_structure_zone,
     &	      vsv_structure_zone,vsh_structure_zone,
     &	      eta_structure_zone,
     &	      qkappa_structure_zone,qmu_structure_zone,
     &	      grid_rho,grid_Ak,grid_Am,grid_Ck,grid_Cm,
     &	      grid_L,grid_N,grid_Fk,grid_Fm,
     &	      grid_kappa,grid_mu,grid_qkappa,grid_qmu,
     &	      idim1_sph,idim2_sph,idim1_tor,idim2_tor )
	  call assign_source
     &	    ( ngrid_r,grid_r,source_r,
     &	      grid_L,grid_Ck,grid_Cm,grid_Fk,grid_Fm,
     &	      grid_qkappa,grid_qmu,grid_mu,
     &	      igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &	      grid_qkappas,grid_qmus,
     &	      idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0 )
c	  call assign_station
c     &	    ( ngrid_r,grid_r,grid_mu,
c     &	      idim_station_sph,idim_station_tor,
c     &	      idim_station_sph0,idim_station_tor0 )
c WENBO
          call assign_station
     &     ( ngrid_r,grid_r,grid_mu,
     &       idim_sph_forcut,idim_tor_forcut,
     &       idim_sph_forcut0,idim_tor_forcut0,
     &       minidim_sph_sta,minidim_tor_sta,
     &       minidim_sph_sta0,minidim_tor_sta0,
     &       idim_ir_sph,idim_ir_tor,idim_ir_sph0,
     &       ir_dep_solid,ir_dep_fluid,ndep_solid,
     &       ndep_fluid,idep_forcut,min_ir_dep)

c **********************************************************************
c computing submatrices for the elastic part of the medium
c **********************************************************************
	  call comp_submatrix
     &	    ( ngrid_r,grid_r,
     &	      grid_rho,grid_kappa,grid_mu,grid_Ak,grid_Am,
     &	      grid_Ck,grid_Cm,grid_L,grid_N,grid_Fk,grid_Fm,
     &	      submatrix_I0,submatrix_I1k,submatrix_I1m,submatrix_I2,
     &	      submatrix_I3k,submatrix_I3m,submatrix_I4,submatrix_I5k,
     &	      submatrix_I5m,submatrix_I6,submatrix_I7 )
	  call comp_submatrix_mod
     &	    ( ngrid_r,grid_r,
     &	      grid_rho,grid_kappa,grid_mu,
     &	      grid_Ak,grid_Am,grid_L,grid_N,grid_Fk,grid_Fm,
     &	      submatrix_I0,submatrix_I1k,submatrix_I3k,
     &	      submatrix_I3m,submatrix_I4,submatrix_I5k,
     &	      submatrix_I5m,submatrix_I6,submatrix_I7,
     &	      submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod )
c **********************************************************************
c **********************************************************************
c **********************************************************************
c computing wavefield for each frequency (for complex omega)
c **********************************************************************
c **********************************************************************
c **********************************************************************
	  do 200 i_frequency=i_frequency_min(i_frequency_band,
     &	                                     n_frequency_band,
     &	                                     n_frequency)
     &	                     +rank,
     &	                     i_frequency_max(i_frequency_band,
     &	                                     n_frequency_band,
     &	                                     n_frequency),
     &	                     size
c	    write(6,*) i_frequency,i_frequency_band,
c     &	               omega(time_series_length,i_frequency,omega_imag)
           print *,'ifreq_start',i_frequency
           lmax_computed=0
	   if ( i_frequency.ne.0 ) then
c	    if ( i_frequency.gt.maxnfreq ) stop 'Too large i_freqnecy.'
            call init_complex_array
     &        ( (maxlmax+1)*5*2*3,coef_c )
            call init_complex_array
     &        ( (maxlmax+1)*5*2*3,coef_dcdr )
	    call init_complex_array
     &	      ( 3*max_nstation,station_disp_solid )
            call init_complex_array
     &        ( 3*max_nstation,station_disp_fluid )
            call init_complex_array
     &        (   max_nstation,epsilon11 )
            call init_complex_array
     &        (   max_nstation,epsilon22 )
            call init_complex_array
     &        (   max_nstation,epsilon33 )
            call init_complex_array
     &        (   max_nstation,epsilon23 )
            call init_complex_array
     &        (   max_nstation,epsilon13 )
            call init_complex_array
     &        (   max_nstation,epsilon12 )
            call init_complex_array
     &        (   max_nstation,pressure )


c
	    do 110 l=lmin(i_frequency_band,n_frequency_band,lmin0 ),
     &	             lmax(i_frequency_band,n_frequency_band,lmax0 )
c **********************************************************************
c computing the displacement and the traction at the boundary
c of the heterogeneous region.
c **********************************************************************
	      if ( l.eq.0 ) then
c WENBO 
               if(source_type.eq.1) then
	        call comp_excitation0
     &	        ( maxngrid_r,
     &	          omega(time_series_length,i_frequency,omega_imag),
     &	          ngrid_r,grid_r,l,source_r,source_mt,
     &	          igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &	          grid_qkappas,grid_qmus,
     &	          submatrix_I0,submatrix_I1k,submatrix_I1m,
     &	          submatrix_I2,submatrix_I3k,submatrix_I3m,
     &	          submatrix_I4,submatrix_I5k,submatrix_I5m,
     &	          submatrix_I6,submatrix_I7,
     &	          submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &	          idim_rs_sph0,idim_rs_tor0,
     &	          whole_vector_sph,whole_vector_tor )
                else
                 call comp_excitation0_single_force
     &           ( maxngrid_r,l,idim_rs_sph0,fr,whole_vector_sph )
                end if
c WENBO
 
	        call comp_wavefield0
     &	        ( maxngrid_r,
     &	          omega(time_series_length,i_frequency,omega_imag),
     &	          submatrix_I0,submatrix_I1k,submatrix_I1m,
     &	          submatrix_I2,submatrix_I3k,submatrix_I3m,
     &	          submatrix_I4,submatrix_I5k,submatrix_I5m,
     &	          submatrix_I6,submatrix_I7,
     &	          submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &	          ngrid_r,grid_r,grid_mu,grid_qkappa,grid_qmu,l,
     &	          idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &	          idim0,init_npos_sph,init_npos_tor,
     &	          idim_rs_sph0,idim_rs_tor0,
     &	          minidim_sph_sta0,minidim_tor_sta0,
     &	          whole_matrix_sph,whole_matrix_tor,
     &	          whole_matrix_dr_sph,whole_matrix_dr_tor,
     &	          whole_vector_sph,whole_vector_tor,work_vector )
	        call check_amp_significance
     &	        ( maxngrid_r,l,nl_check_amp,
     &            amp_max,i_significance,idim0,
     &	          init_npos_sph,init_npos_tor,i_frequency,
     &	          idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &	          grid_mu,idim_sph_forcut0,idim_tor_forcut0,
     &	          whole_vector_sph,whole_vector_tor,work_vector,
     &            min_ir_dep)
c	        call comp_displacement_station0
c     &	        ( maxngrid_r,maxlmax,
c     &	          whole_vector_tor,whole_vector_sph,
c     &	          l,n_station,station_theta,station_phi,
c     &	          idim_station_sph0,idim_station_tor0,
c     &	          vecsph_sph1,vecsph_sph2,vecsph_tor,
c     &	          station_displacement )

                call comp_disp_strain_all0(station_disp_solid,
     &            station_disp_fluid,
     &            epsilon11,epsilon22,epsilon33,epsilon12,
     &            epsilon13,epsilon23,pressure,
     &            nl_check_amp,i_significance,
     &            maxlmax,coef_c,coef_dcdr,maxngrid_r,max_ndep,
     &            r_sta_solid,r_sta_fluid,
     &            ir_dep_solid,ir_dep_fluid,idim_ir_sph0,
     &            ir_for_stress,ir_for_pressure,
     &            ir_top_stress,ir_bot_stress,
     &            ir_top_pressure,ir_bot_pressure,
     &            whole_vector_sph,l,ir_CMB,ir_ICB,grid_rho,
     &            n_station_solid,n_station_fluid,
     &            ndist_solid,ndist_fluid,dist_solid,dist_fluid,
     &            omega(time_series_length,i_frequency,omega_imag),
     &            grid_r,ngrid_r,
     &            depth_solid,depth_fluid,
     &            ndep_solid,ndep_fluid,top_fluid,bot_fluid,
     &            source_type,i_frequency)
                  lmax_computed=0
	      else
c WENBO
c	        if ( i_significance.eq.1 ) then
c               i_significance=1: not convergent, continue calculation.
c               i_significance=0: convergent, 
c                     but end calculation in next round(i_significance=-1). 
c               i_significance=-1: calculation finished
                if ( i_significance.ge.0 ) then

c WENBO
                 if(l.eq.lmax(i_frequency_band,n_frequency_band,lmax0)
     &              .and.i_significance.ge.0) then
                    write(6,*) "WARNING:Not reach convergence!!!!!!!"
                 end if

                 if (source_type.eq.1) then
	          call comp_excitation
     &	          ( maxngrid_r,
     &	            omega(time_series_length,i_frequency,omega_imag),
     &	            ngrid_r,grid_r,l,source_r,source_mt,
     &	            igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &	            grid_qkappas,grid_qmus,
     &	            submatrix_I0,submatrix_I1k,submatrix_I1m,
     &	            submatrix_I2,submatrix_I3k,submatrix_I3m,
     &	            submatrix_I4,submatrix_I5k,submatrix_I5m,
     &	            submatrix_I6,submatrix_I7,
     &	            submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &	            idim_rs_sph,idim_rs_tor,
     &	            whole_vector_sph,whole_vector_tor )
                  else
                   call comp_excitation_single_force
     &              (maxngrid_r,l,idim_rs_sph,idim_rs_tor,fr,ftheta,
     &               fphi,whole_vector_sph,whole_vector_tor)
                  end if

	        call comp_wavefield
     &	        ( maxngrid_r,
     &	          omega(time_series_length,i_frequency,omega_imag),
     &	          submatrix_I0,submatrix_I1k,submatrix_I1m,
     &	          submatrix_I2,submatrix_I3k,submatrix_I3m,
     &	          submatrix_I4,submatrix_I5k,submatrix_I5m,
     &	          submatrix_I6,submatrix_I7,
     &	          submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &	          ngrid_r,grid_r,grid_mu,grid_qkappa,grid_qmu,l,
     &	          idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &	          idim0,init_npos_sph,init_npos_tor,
     &	          idim_rs_sph,idim_rs_tor,
     &	          minidim_sph_sta,minidim_tor_sta,
     &	          whole_matrix_sph,whole_matrix_tor,
     &	          whole_matrix_dr_sph,whole_matrix_dr_tor,
     &	          whole_vector_sph,whole_vector_tor,work_vector,
     &            source_type )

c	        call comp_displacement_station
c     &	        ( maxngrid_r,maxlmax,
c     &	          whole_vector_tor,whole_vector_sph,
c     &	          l,n_station,station_theta,station_phi,
c     &	          idim_station_sph,idim_station_tor,
c     &	          vecsph_sph1,vecsph_sph2,vecsph_tor,
c     &	          station_displacement(1,1,i_frequency) )

                 call comp_disp_strain_all(station_disp_solid,
     &            station_disp_fluid,
     &            epsilon11,epsilon22,epsilon33,epsilon12,
     &            epsilon13,epsilon23,pressure,
     &            maxlmax,coef_c,coef_dcdr,maxngrid_r,max_ndep,
     &            r_sta_solid,r_sta_fluid,
     &            ir_dep_solid,ir_dep_fluid,ir_for_stress,
     &            ir_for_pressure,ir_top_stress,ir_bot_stress,
     &            ir_top_pressure,ir_bot_pressure,
     &            idim_ir_sph,idim_ir_tor,
     &            whole_vector_tor,whole_vector_sph,
     &            nl_check_amp,i_significance,
     &            l,ir_CMB,ir_ICB,grid_rho,
     &            n_station_solid,n_station_fluid,
     &            ndist_solid,ndist_fluid,dist_solid,dist_fluid,rank,
     &            omega(time_series_length,i_frequency,omega_imag),
     &            grid_r,ngrid_r,depth_solid,depth_fluid,
     &            ndep_solid,ndep_fluid,top_fluid,bot_fluid,
     &            source_type,i_frequency,pi)


!WENBO
                call check_amp_significance
     &          ( maxngrid_r,l,nl_check_amp,
     &            amp_max,i_significance,idim0,
     &            init_npos_sph,init_npos_tor,i_frequency,
     &            idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &            grid_mu,idim_sph_forcut,idim_tor_forcut,
     &            whole_vector_sph,whole_vector_tor,work_vector,
     &            min_ir_dep)

                  lmax_computed=l
                 if(mod(l,2000).eq.0) print *,'l=',l,'ifreq=',
     &              i_frequency,'-done'
	        endif
	      endif
  110	    continue
	   endif
c **********************************************************************
c writing the result to the spectrum files
c **********************************************************************
        call write_disp_velo_stress(i_frequency,
     &     time_series_length,omega_imag,
     &     ngrid_r,n_station_solid,n_station_fluid,
     &     r_sta_solid,r_sta_fluid,
     &     ir_dep_solid,max_ndep,maxn_structure_zone,
     &     izone_idep_solid,izone_idep_fluid,
     &     n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &     rho_structure_zone,vpv_structure_zone,vph_structure_zone,
     &     vsv_structure_zone,vsh_structure_zone,eta_structure_zone,
     &     save_velo,ndep_solid,ndist_solid,
     &     max_nstation,station_disp_solid,
     &     station_disp_fluid,pressure,
     &     epsilon11,epsilon22,epsilon33,
     &     epsilon12,epsilon13,epsilon23,
     &     source_type,top_fluid,bot_fluid,lmax_computed,
     &     maxlmax,coef_c,coef_dcdr,Atop,Ctop,Ftop,Ltop,Ntop,
     &     Abot,Cbot,Fbot,Lbot,Nbot)

!ndep_fluid,ntheta_fluid,

c	   call write_sac_file
c     &	      ( i_frequency,n_station,sac_file,
c     &	        station_displacement )
           print *,'ifreq=',i_frequency,' done'
  200	  continue
c
 1000	continue
c
c	  nsnd = 2 * 3 * max_nstation
c	  if ( rank.ne.0 ) then
c	    call msndi
c	    do 1030 i_frequency_band=1,n_frequency_band
c	    do 1020 i_frequency=i_frequency_min(i_frequency_band,
c     &	                                       n_frequency_band,
c     &	                                       n_frequency)
c     &	                       +rank,
c     &	                       i_frequency_max(i_frequency_band,
c     &	                                       n_frequency_band,
c     &	                                       n_frequency),
c     &	                       size
c	      call mputd( station_displacement(1,1,i_frequency),
c     &	                  nsnd,info )
 1020	    continue
 1030	    continue
c	    call msend( 0,20,info )
c	  else
c	    do 1060 irank=1,size-1
c	      call mrecv( 20,irank,info )
c	      do 1050 i_frequency_band=1,n_frequency_band
c	      do 1040 i_frequency=i_frequency_min(i_frequency_band,
c     &	                                       n_frequency_band,
c     &	                                       n_frequency)
c     &	                       +irank,
c     &	                       i_frequency_max(i_frequency_band,
c     &	                                       n_frequency_band,
c     &	                                       n_frequency),
c     &	                       size
c	        call mgetd( station_displacement(1,1,i_frequency),
c     &	                    nsnd,info )
 1040	      continue
 1050	      continue
 1060	    continue
c	  endif
 1100	continue
c	if ( rank.eq.0 ) then
c	  call convspc
c     &	    ( max_nstation,maxnfreq,n_station,
c     &	   time_series_length,n_frequency,omega_imag,
c     &	   station_displacement,
c     &	   source_depth,source_lat,source_lon,station_lat,station_lon,
c     &	   work_spc,work_time,sac_file )
c	endif

        call mpi_reduce(lmax_allfreq,lmax_computed,1,mpi_integer,
     &       MPI_MAX,0,mpi_comm_world,ierr )


c          lmax_allfreq=lmax_computed


        if(rank.eq.0) call  save_par(n_frequency,rank,
     &         ndist_solid,dist_solid,ndist_fluid,
     &         dist_fluid,max_ntheta,
     &         ndep_solid,ndep_fluid,source_type,
     &         time_series_length,omega_imag,lmax_allfreq,
     &         top_fluid,bot_fluid,ir_dep_fluid,ir_dep_solid,
     &         max_ndep,maxngrid_r,grid_r,grid_rho,
     &         Atop,Ctop,Ftop,Ltop,Ntop,
     &         Abot,Cbot,Fbot,Lbot,Nbot)
        if(rank.eq.0)print *,'finally done!'


c
	call mpi_finalize( info )
c
	end
