cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine param_input
     &	  ( maxngrid_r,maxlmax,max_nstation,maxn_structure_zone,
     &	    time_series_length,n_frequency,omega_imag,
     &	    max_ndep,max_ntheta,ngrid_r,lmin,lmax,
     &      n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &      rho_structure_zone,vpv_structure_zone,vph_structure_zone,
     &      vsv_structure_zone,vsh_structure_zone,eta_structure_zone,
     &      qkappa_structure_zone,qmu_structure_zone,
     &      fluid_thiszone,r_CMB,r_ICB,
     &      source_r,source_type,source_mt,fr,ftheta,fphi,
     &      source_depth,source_lat,source_lon,
     &      depth_solid,depth_fluid,ndist_ela,
     &      ndist_acou,dist_elastic,dist_acoustic,
     &      ndep_solid,ndep_fluid,depth_for_stress,depth_for_pressure,
     &      izone_idep_solid,izone_idep_fluid,
     &      depth_tolerence,top_fluid,
     &      bot_fluid,idep_forcut,
     &      save_velo,n_station_solid,n_station_fluid)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c inputting parameters
c    required subroutines: error_handling,distaz
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
        implicit none
        integer maxngrid_r,maxlmax,max_nstation,maxn_structure_zone
        integer max_ndep,max_ntheta
	integer n_frequency,ngrid_r,lmin,lmax
	integer n_structure_zone
        integer n_station_solid,n_station_fluid
	real*8 time_series_length,omega_imag
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
        integer fluid_thiszone(maxn_structure_zone)
        real*8 r_CMB,r_ICB
	real*8 source_r,source_mt(3,3)
c WENBO
        real*8 fr,ftheta,fphi
        integer source_type
c WENBO
	real*8 source_depth,source_lat,source_lon
        real*8 Ddepth
c        real*8 min_dep,min_dep_solid,min_dep_fluid,ddep
c        integer ndep,ndep_solid,ndep_fluid,idep_forcut
        integer ndep_solid,ndep_fluid,idep_forcut
        integer top_fluid,bot_fluid
        real*8 depth_solid(max_ndep),depth_fluid(max_ndep)
        integer izone_idep_solid(max_ndep)
        integer izone_idep_fluid(max_ndep)
        real*8  depth_for_stress(3,max_ndep)
        real*8  depth_for_pressure(3,max_ndep)
        real*8  depth_tolerence,ddepth_for_interpo
        integer ndist_ela,ndist_acou
        real*8 dist_elastic(max_ntheta),dist_acoustic(max_ntheta)
        integer save_velo
c	real*8 station_lat(max_nstation),station_lon(max_nstation)
c	real*8 station_theta(max_nstation),station_phi(max_nstation)
c	character*80 sac_file(max_nstation)
c other variables
        integer ndist_ela_tmp,ndist_acou_tmp
        real*8 dist_elastic_tmp,dist_acoustic_tmp,depth_tole_tmp
	integer i,igll,inode,idep,idep_recount,idist,istat,itheta
        integer nexp,nerr,nwgll,izone_thisdep
	real dist,az,baz,xdeg,min_dep,max_dep
	character*80 dummy,tmpfile,tmpfile1,dep_solid_file,
     &                dep_fluid_file,dist_solid_file,dist_fluid_file
        real*8 r_up,r_down,r_freesurf
c constants
	real*8 pi
	parameter ( pi=3.1415926535897932d0 )
c
	data tmpfile  / 'work_r1' /
        data tmpfile1 / 'work_tmp' /
c
c **********************************************************************
c reading input file from the standard input and writing out to the
c temporary file
c **********************************************************************
c opening a temporary file
	open( unit=11, file=tmpfile1, status='unknown' )
c writing to the temporary file
  100	continue
	  read(5,110) dummy
  110	  format(a80)
	  if ( dummy(1:1).eq.'c' ) goto 100
	  if ( dummy(1:3).eq.'end' ) goto 120
	  write(11,110) dummy
	  goto 100
  120	continue
c closing the temporary file
	close(11)
c
c **********************************************************************
c reading the parameter from the temporary file
c **********************************************************************
c opening the temporary file
	open( unit=11, file=tmpfile1, status='unknown' )
c reading the parameters
c ---- parameters for time series ---
	read(11,*) time_series_length,n_frequency
	read(11,*) omega_imag
c ---- parameters for numerical grids ----
	read(11,*) ngrid_r,lmin,lmax
	if ( ngrid_r.gt.maxngrid_r )
     &	  call error_handling(1)
c        if ( lmax.gt.maxlmax ) print *,lmax,maxlmax
	if ( lmax.gt.maxlmax )
     &	  call error_handling(3)
c ---- parameters for structure ---
	read(11,*) n_structure_zone
	if ( n_structure_zone.gt.maxn_structure_zone )
     &	  call error_handling(2)
        r_CMB=0.d0
        r_ICB=0.d0
	do 130 i=1,n_structure_zone
	  read(11,*) rmin_structure_zone(i),rmax_structure_zone(i),
     &	             rho_structure_zone(1,i),rho_structure_zone(2,i),
     &	             rho_structure_zone(3,i),rho_structure_zone(4,i)
	  read(11,*) vpv_structure_zone(1,i),vpv_structure_zone(2,i),
     &	             vpv_structure_zone(3,i),vpv_structure_zone(4,i)
	  read(11,*) vph_structure_zone(1,i),vph_structure_zone(2,i),
     &	             vph_structure_zone(3,i),vph_structure_zone(4,i)
	  read(11,*) vsv_structure_zone(1,i),vsv_structure_zone(2,i),
     &	             vsv_structure_zone(3,i),vsv_structure_zone(4,i)
	  read(11,*) vsh_structure_zone(1,i),vsh_structure_zone(2,i),
     &	             vsh_structure_zone(3,i),vsh_structure_zone(4,i)
	  read(11,*) eta_structure_zone(1,i),eta_structure_zone(2,i),
     &	             eta_structure_zone(3,i),eta_structure_zone(4,i),
     &	             qmu_structure_zone(i),qkappa_structure_zone(i)

          if ( ( ( vsv_structure_zone(1,i).eq.0.d0 ).and.
     &           ( vsv_structure_zone(2,i).eq.0.d0 ).and.
     &           ( vsv_structure_zone(3,i).eq.0.d0 ).and.
     &           ( vsv_structure_zone(4,i).eq.0.d0 )      ).or.
     &         ( ( vsh_structure_zone(1,i).eq.0.d0 ).and.
     &           ( vsh_structure_zone(2,i).eq.0.d0 ).and.
     &           ( vsh_structure_zone(3,i).eq.0.d0 ).and.
     &           ( vsh_structure_zone(4,i).eq.0.d0 )    ) ) then
                fluid_thiszone(i)=1
          else
                fluid_thiszone(i)=0
          end if

          if(i.ge.2.and.i.lt.n_structure_zone) then
              if(dabs(vsv_structure_zone(1,i)).lt.1.e-7.and.
     &           dabs(vsv_structure_zone(1,i-1)).gt.1.e-7) then
                  r_ICB=rmin_structure_zone(i)
              end if
             if(dabs(vsv_structure_zone(1,i)).gt.1.e-7.and.
     &          dabs(vsv_structure_zone(1,i-1)).lt.1.e-7) then
                  r_CMB=rmin_structure_zone(i)
             end if
               
          end if
  130	continue

        r_freesurf=rmax_structure_zone(n_structure_zone)
        if(r_CMB.lt.1.e-7.or.r_ICB.lt.1.e-7) then
c               print *,'r_CMB',r_CMB,r_ICB
               stop 'Error in finding r_CMB or r_ICB'
        end if

c ---- parameters for a source ---
	read(11,*) source_depth,source_lat,source_lon,source_type
c WENBO
        if(source_type.eq.1) then
	   read(11,*) nexp,
     &	           source_mt(1,1),source_mt(2,2),source_mt(3,3),
     &	           source_mt(1,2),source_mt(1,3),source_mt(2,3)
        else if(source_type.eq.2) then
           read(11,*) nexp,
     &             fr,ftheta,fphi
        else
           stop 'Error of source_type'
        end if
c WENBO
	source_r = rmax_structure_zone(n_structure_zone)
     &	            - source_depth
	source_mt(1,1) = source_mt(1,1) * ( 10.d0**(nexp-25) )
	source_mt(2,2) = source_mt(2,2) * ( 10.d0**(nexp-25) )
	source_mt(3,3) = source_mt(3,3) * ( 10.d0**(nexp-25) )
	source_mt(1,2) = source_mt(1,2) * ( 10.d0**(nexp-25) )
	source_mt(1,3) = source_mt(1,3) * ( 10.d0**(nexp-25) )
	source_mt(2,3) = source_mt(2,3) * ( 10.d0**(nexp-25) )
	source_mt(2,1) = source_mt(1,2)
	source_mt(3,1) = source_mt(1,3)
	source_mt(3,2) = source_mt(2,3)
        fr=fr*( 10.d0**(nexp-20) )
        ftheta=ftheta*( 10.d0**(nexp-20) )
        fphi=fphi*( 10.d0**(nexp-20) )

c --- parameters for stations ---
        read(11,*)dep_solid_file
        read(11,*)dist_solid_file
        read(11,*)dep_fluid_file
        read(11,*)dist_fluid_file


        read(11,*)save_velo
        close(11)

!        call distaz(0.0,0.0,0.0,90.0,1,dist,az,baz,xdeg,nerr)

        open( unit=12, file=dep_solid_file, status='unknown' )
!The pamarter used to find and/or add nodes for stress calculation.
        read(12,*) ddepth_for_interpo,ndep_solid

!The parameter used to check whether adding another node or two nodes collopsed to one.
!Therefore depth_tole_tmp have to be larger than ddepth_for_interpo
        read(12,*) depth_tole_tmp
        read(12,*)nwgll
        if(nwgll.lt.2.or.nwgll.gt.10)
     &         stop 'nwgll should be between 2 and 10'

        depth_tolerence=depth_tole_tmp
        if(ddepth_for_interpo.lt.depth_tolerence) 
     &      stop 'Error,ddepth_for_interpo.lt.depth_tolerence'
        if(ndep_solid.gt.max_ndep) stop 'Error,ndep_solid>max_ndep'
!        if(ndep_solid.le.2) stop 'ndep_solid should be larger than 2'
        idep_recount=0
        do idep=1,ndep_solid
           read(12,*)depth_solid(idep),izone_idep_solid(idep)
        end do
        
        do idep=1,ndep_solid
           idep_recount=idep_recount+1
           igll=mod(idep_recount,nwgll-1)
           if(igll.eq.0) igll=nwgll-1

           if(idep.eq.1) then
             depth_for_stress(3,idep)=depth_solid(idep) !up
             depth_for_stress(2,idep)=depth_solid(idep)+
     &           ddepth_for_interpo  !mid
             depth_for_stress(1,idep)=depth_solid(idep)+
     &           2.0*ddepth_for_interpo  !down
           else if(izone_idep_solid(idep).ne.
c discontinuity
     &             izone_idep_solid(idep-1)) then
!update the depths (between idep and idep-1) setten in the previous step  
             depth_for_stress(3,idep)=depth_solid(idep) !up
             depth_for_stress(2,idep)=depth_solid(idep)+
     &               ddepth_for_interpo    !mid
             depth_for_stress(1,idep)=depth_solid(idep)+
     &               2.0*ddepth_for_interpo    !down

!depth between idep and idep+1
             depth_for_stress(3,idep-1)=depth_solid(idep-1) -
     &               2.0*ddepth_for_interpo    !up

             depth_for_stress(2,idep-1)=depth_solid(idep-1) -
     &               ddepth_for_interpo    !mid
             depth_for_stress(1,idep-1)=depth_solid(idep-1) !down 

             idep_recount=idep_recount-1
           else if(idep.eq.ndep_solid) then
             depth_for_stress(3,idep)=depth_solid(idep)  -
     &               2.0*ddepth_for_interpo    !up
             depth_for_stress(2,idep)=depth_solid(idep) -
     &               ddepth_for_interpo    !mid
             depth_for_stress(1,idep)=depth_solid(idep) !down
!to be fixed, it only works for HIGH_RESOLUTION (all the gll points
!involved).
!             if(igll.ne.1) stop 'Error of the depth table!'
           else 
              depth_for_stress(3,idep)= depth_solid(idep) -
     &               ddepth_for_interpo    !up
              depth_for_stress(2,idep)= depth_solid(idep) !mid
              depth_for_stress(1,idep)= depth_solid(idep) +
     &               ddepth_for_interpo    !down
           end if
!check the validaty
           izone_thisdep=izone_idep_solid(idep)

              if( r_freesurf-depth_solid(idep)-depth_tolerence.gt.
     &            rmax_structure_zone(izone_thisdep).or.
     &            r_freesurf- depth_solid(idep)+depth_tolerence.lt.
     &            rmin_structure_zone(izone_thisdep) )
     &             stop 'Depth and izone are not consistent'

           if(idep.gt.2) then

              if(depth_solid(idep).lt.depth_solid(idep-1).and.
     &           izone_idep_solid(idep).eq.izone_idep_solid(idep-1)) 
     &          stop 'depth_solid table should be ascending'

              if(dabs(depth_solid(idep)-depth_solid(idep-1)).lt.
     &             depth_tolerence/2.d0 .and. 
     &           izone_idep_solid(idep).eq.izone_idep_solid(idep-1))  
     &           stop 'Depth difference is too small or 
     &                depth_tolerence is too large! That makes
     &                finding ir_dep_stress impossible!'
           end if
        end do
        close(12)




        open( unit=13, file=dist_solid_file, status='unknown' )
        istat=0
        do idep=1,ndep_solid
           read(13,*) ndist_ela_tmp
           if(idep.eq.2) ndist_ela=ndist_ela_tmp
           itheta=0
           do idist=1,ndist_ela_tmp
             itheta=itheta+1
             read(13,*) dist_elastic_tmp
             if(idep.eq.2) dist_elastic(itheta)=dist_elastic_tmp*pi/180.0
           end do
        end do
        close(13)


        open( unit=14, file=dep_fluid_file, status='unknown' )
        read(14,*) ddepth_for_interpo,ndep_fluid
        read(14,*) depth_tole_tmp
        depth_tolerence=depth_tole_tmp
        if(depth_tole_tmp.ne.depth_tolerence) 
     &       stop  'depth_tolerence for fluid and solid are different'
        if(ndep_fluid.gt.max_ndep) stop 'Error,ndep_fluid>max_ndep'

        read(14,*)nwgll
        if(nwgll.lt.2.or.nwgll.gt.10)
     &         stop 'nwgll should be between 2 and 10'
        idep_recount=0
        do idep=1,ndep_fluid
           read(14,*)depth_fluid(idep),izone_idep_fluid(idep)
           idep_recount=idep_recount+1

           igll=mod(idep_recount,nwgll-1)
!           if(igll.eq.1) igll=nwgll
           if(igll.eq.0) igll=nwgll-1

           if(idep.eq.1) then
             depth_for_pressure(3,idep)=depth_fluid(idep) !up
             depth_for_pressure(2,idep)=depth_fluid(idep)+
     &           ddepth_for_interpo  !mid
             depth_for_pressure(1,idep)=depth_fluid(idep)+
     &           2.0*ddepth_for_interpo  !down
           else if(izone_idep_fluid(idep).ne.
     &             izone_idep_fluid(idep-1)) then
             depth_for_pressure(3,idep)=depth_fluid(idep) !up
             depth_for_pressure(2,idep)=depth_fluid(idep)+
     &               ddepth_for_interpo    !mid
             depth_for_pressure(1,idep)=depth_fluid(idep)+
     &               2.0*ddepth_for_interpo    !down
! the upper idepth
             depth_for_pressure(1,idep-1)=depth_fluid(idep)!down
             depth_for_pressure(2,idep-1)=depth_fluid(idep)-
     &               ddepth_for_interpo    !mid
             depth_for_pressure(3,idep-1)=depth_fluid(idep)-
     &               2.0*ddepth_for_interpo    !up

             idep_recount=idep_recount-1
           else if(idep.eq.ndep_fluid) then
             depth_for_pressure(3,idep)=depth_fluid(idep)  -
     &               2.0*ddepth_for_interpo    !up
             depth_for_pressure(2,idep)=depth_fluid(idep)  -
     &               ddepth_for_interpo    !mid
             depth_for_pressure(1,idep)=depth_fluid(idep)  !down
            if(igll.ne.1) stop 'Error of the depth table (fluid)!'
           else
              depth_for_pressure(3,idep)=depth_fluid(idep)  -
     &               ddepth_for_interpo    !up
              depth_for_pressure(2,idep)=depth_fluid(idep) !mid
              depth_for_pressure(1,idep)=depth_fluid(idep) +
     &               ddepth_for_interpo    !down

           end if
!check the validaty
           izone_thisdep=izone_idep_fluid(idep)
              if( r_freesurf-depth_fluid(idep)-depth_tolerence.gt.
     &            rmax_structure_zone(izone_thisdep).or.
     &            r_freesurf- depth_fluid(idep)+depth_tolerence.lt.
     &            rmin_structure_zone(izone_thisdep) )
     &            stop 'Depth(fluid) and izone are not consistent'

           if(idep.gt.2) then
              if(depth_fluid(idep).lt.depth_fluid(idep-1).and.
     &           izone_idep_fluid(idep).eq.
     &           izone_idep_fluid(idep-1))
     &          stop 'depth_fluid table should be ascending'

              if(dabs(depth_fluid(idep)-depth_fluid(idep-1)).lt.
     &             depth_tolerence/2.d0 .and.
     &           izone_idep_fluid(idep).eq.
     &           izone_idep_fluid(idep-1) )
     &           stop 'Depth difference is too small or 
     &                depth_tolerence is too large! That makes
     &                finding ir_dep_pressure impossible!'
           end if

        end do
        close(14)
        
        if(ndep_fluid.lt.2.and.ndep_solid.lt.2) 
     &    stop 'Either ndep_fluid or ndep_solid should be larger than 2'

        open( unit=15, file=dist_fluid_file, status='unknown' )
        istat=0
        do idep=1,ndep_fluid
           read(15,*) ndist_acou_tmp
           if(idep.eq.2) then
               ndist_acou=ndist_acou_tmp
               if(ndist_acou.ne.ndist_ela.and.ndep_solid.gt.0) 
     &           stop 'ndist_acou should be equal to ndist_ela'
           end if
           itheta=0
           do idist=1,ndist_acou_tmp
             itheta=itheta+1
             read(15,*) dist_acoustic_tmp
             if(idep.eq.2) then
               dist_acoustic(itheta)=dist_acoustic_tmp*pi/180.0
               if(dabs(dist_acoustic(itheta)-dist_elastic(itheta)).gt.
     &            1.e-10.and.ndep_solid.gt.0) then
                stop 'distance table of elastic and acoustic 
     &                 should be  exactly the same'
               end if
             end if
           end do
        end do
        close(15)


        if(ndep_fluid.le.0.and.ndep_solid.le.0) 
     &      stop 'both ndep_solid and ndep_fluid <0'

        min_dep=1.e10
        max_dep=-1.e10
        Ddepth=1.e12
        idep_forcut=-1

        if(ndep_solid.gt.0) then
          min_dep=depth_solid(1)
          max_dep=depth_solid(ndep_solid)
          do idep=1,ndep_solid
            if(idep.eq.1) then
                 idep_forcut=1
                 Ddepth=dabs(depth_solid(1)-source_depth)
            else if(dabs(depth_solid(idep)-source_depth).lt.Ddepth.and.
     &          Ddepth.ge.25.0) then
                 idep_forcut=idep
                 Ddepth=dabs(depth_solid(idep)-source_depth)
            else if(dabs(depth_solid(idep)-source_depth).gt.Ddepth.and.
     &          Ddepth.lt.25.0) then
                 idep_forcut=idep
                 Ddepth=dabs(depth_solid(idep)-source_depth)
            end if
          end do
        end if
       
        if(ndep_fluid.gt.0) then
          if(min_dep.gt.depth_fluid(1)) then
              min_dep=depth_fluid(1)
              top_fluid=1
          else
              top_fluid=0
          end if
          if(max_dep.lt.depth_fluid(ndep_fluid)) then
              max_dep=depth_fluid(ndep_fluid)
              bot_fluid=1
          else
              bot_fluid=0
          end if
          do idep=1,ndep_fluid
             if(idep.eq.1.and.ndep_solid.eq.0) then
                 idep_forcut=1
                 Ddepth=dabs(depth_fluid(1)-source_depth)
             else if(dabs(depth_fluid(idep)-source_depth).lt.Ddepth.and.
     &          Ddepth.ge.25.0) then
                 idep_forcut=ndep_solid+idep
                 Ddepth=dabs(depth_fluid(idep)-source_depth)
             else if(dabs(depth_fluid(idep)-source_depth).gt.Ddepth.and.
     &          Ddepth.lt.25.0) then
                 idep_forcut=ndep_solid+idep
                 Ddepth=dabs(depth_fluid(idep)-source_depth)
             end if
          end do
        else
          top_fluid=0
          bot_fluid=0
        end if

c     If the target depth is close to source depth, the convergence
c     would be very slow (Refer to Fig.7 in Takeuchi and Geller, GJI 2006).
c     Thus, we set threshold depth at 20km and using average to reduce
c     the possible undulation for points near source depth.
        if(Ddepth.lt.25.0) then
          print *,'WARNING!!!!!!!! The depth_cut_off is close to 
     &      depth_source (dif_depth<5km). Thus, the convergence 
     &      rate is slow (too many modes (l) get involved)
     &      and it wastes compuation resources!!!!'
        end if
        r_up=rmax_structure_zone(n_structure_zone)-min_dep
        r_down=rmax_structure_zone(n_structure_zone)-max_dep
        if(r_up>r_CMB.and.r_down<r_ICB) then
             stop 'Error, the depth range from dep_min to 
     &            dep_max crosses the whole outer core,
     &            that is not allowed in present code!'
        end if


        n_station_solid=ndep_solid*ndist_ela
        n_station_fluid=ndep_fluid*ndist_acou



!        call distaz(0.0,0.0,0.0,90.0,1,dist,az,baz,xdeg,nerr)
!        print *,'az',az

        if ( n_station_solid.gt.max_nstation.or.
     &        n_station_fluid.gt.max_nstation ) call error_handling(5)


	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine grid_generation
     &	  ( maxngrid_r,ngrid_r0,
     &	    n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &	    vsv_structure_zone,vsh_structure_zone,
     &      fluid_thiszone,r_CMB,r_ICB,ir_CMB,ir_ICB,
     &	    grid_r,source_r,ngrid_r,depth_solid,depth_fluid,
     &      ndep_solid,ndep_fluid,max_ndep,depth_for_stress,
     &      depth_for_pressure,ir_for_stress,ir_for_pressure,
     &      ir_dep_solid,ir_dep_fluid,r_sta_solid,r_sta_fluid,
     &      depth_tolerence,izone_idep_solid,izone_idep_fluid,
     &      r_top_stress,r_bot_stress,
     &      r_top_pressure,r_bot_pressure,
     &      ir_top_stress,ir_bot_stress,ir_top_pressure,
     &      ir_bot_pressure,top_fluid,bot_fluid)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the number and the location of grid points.
c    required subroutines: error_handling
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
        implicit none
	integer maxngrid_r,ngrid_r0,n_structure_zone,ngrid_r
        integer fluid_thiszone
	real*8 rmin_structure_zone(*),rmax_structure_zone(*)
	real*8 vsv_structure_zone(4,*)
	real*8 vsh_structure_zone(4,*)
        real*8 r_CMB,r_ICB
	real*8 grid_r(*),source_r
c WENBO
c        real*8 min_dep_solid,min_dep_fluid,ddep
        integer izone_idep_solid(*),izone_idep_fluid(*)
        real*8 depth_solid(*),depth_fluid(*)
        real*8 depth_for_stress(3,*),depth_for_pressure(3,*)
        real*8 r_top_stress(3),r_bot_stress(3)
        real*8 r_top_pressure(3),r_bot_pressure(3) 
        integer ir_dep_solid(*),ir_dep_fluid(*)
        integer ir_for_stress(3,*),ir_for_pressure(3,*)
        integer ir_top_stress(3),ir_bot_stress(3)
        integer ir_top_pressure(3),ir_bot_pressure(3)
        real*8 r_sta_solid(*),r_sta_fluid(*)
        real*8 depth_tolerence
        integer max_ndep,ndep_solid,ndep_fluid
        integer ir_CMB,ir_ICB
        integer top_fluid,bot_fluid
c WENBO
c other variables
	integer ngrid,m_structure_zone,nnodes_more
	integer i_zone,itmp,i,ir
	real*8 rmin,rmax,rh
c WENBO
        real*8 radius,depth
c        real*8 max2_dr,max1_dr
        logical fluid
c
c computing the layer number of the top solid layer
	m_structure_zone = 0
	do 100 i_zone=1,n_structure_zone
	  if ( ( ( vsv_structure_zone(1,i_zone).eq.0.d0 ).and.
     &	         ( vsv_structure_zone(2,i_zone).eq.0.d0 ).and.
     &	         ( vsv_structure_zone(3,i_zone).eq.0.d0 ).and.
     &	         ( vsv_structure_zone(4,i_zone).eq.0.d0 )      ).or.
     &	       ( ( vsh_structure_zone(1,i_zone).eq.0.d0 ).and.
     &	         ( vsh_structure_zone(2,i_zone).eq.0.d0 ).and.
     &	         ( vsh_structure_zone(3,i_zone).eq.0.d0 ).and.
     &	         ( vsh_structure_zone(4,i_zone).eq.0.d0 )      ) ) then
	    continue
	  else
	    m_structure_zone = i_zone
	  endif
  100	continue
c
c computing the number and the location of the grid points
	rmin = rmin_structure_zone(1)
	rmax = rmax_structure_zone(n_structure_zone)
	grid_r(1) = rmin
	itmp = 1
	do 210 i_zone=1,n_structure_zone
          if (fluid_thiszone.eq.1) then
                fluid=.true.
          else
                fluid=.false.
          end if
	  rh = rmax_structure_zone(i_zone) - rmin_structure_zone(i_zone)
	  ngrid = int( dble(ngrid_r0) * rh / ( rmax - rmin ) ) + 1
	  do 200 i=1,ngrid
	    itmp = itmp + 1
	    if ( itmp.gt.maxngrid_r ) call error_handling(11)
	    grid_r(itmp) = rmin_structure_zone(i_zone)
     &	                            + rh * dble(i) / dble( ngrid )
c ------ putting an additional node on the source location
            if ( itmp.gt.2 ) then
              if ( ( grid_r(itmp-1).lt.source_r ).and.
     &             ( source_r.lt.grid_r(itmp) ) ) then
                if ( itmp+1.gt.maxngrid_r ) call error_handling(11)
                grid_r(itmp+1) = grid_r(itmp)
                grid_r(itmp) = source_r
                itmp = itmp + 1
                if ( ( i_zone.eq.m_structure_zone ).and.
     &               ( i.eq.ngrid )                      ) then
                  if ( itmp+1.gt.maxngrid_r ) call error_handling(11)
                  grid_r(itmp+1) = grid_r(itmp)
                  grid_r(itmp) = ( grid_r(itmp) + source_r ) / 2.d0
                  itmp = itmp + 1
                endif
              else
                if ( ( i_zone.eq.m_structure_zone ).and.
     &               ( i.eq.ngrid ).and.
     &               ( source_r.eq.grid_r(itmp) )        ) then
                  if ( itmp+2.gt.maxngrid_r ) call error_handling(11)
                  source_r = grid_r(itmp-1) * 0.001d0
     &                         + grid_r(itmp) * 0.999d0
                  grid_r(itmp+2) = grid_r(itmp)
                  grid_r(itmp+1) = ( grid_r(itmp) + source_r ) / 2.d0
                  grid_r(itmp) = source_r
                  itmp = itmp + 2
                endif
              endif
            endif
c ------
  200	  continue
  210	continue
c recouting the total number of grid points
        ngrid_r = itmp

        call addnodes_besides_topbot(top_fluid,bot_fluid,
     &           n_structure_zone,rmax_structure_zone,
     &           rmin_structure_zone,izone_idep_solid,
     &           izone_idep_fluid,depth_tolerence,
     &           depth_solid,depth_fluid,ndep_solid,
     &           ndep_fluid,depth_for_stress,depth_for_pressure,
     &           ir_for_stress,ir_for_pressure,r_top_stress,
     &           r_bot_stress,r_top_pressure,r_bot_pressure,
     &           ir_top_stress,ir_bot_stress,
     &           ir_top_pressure,ir_bot_pressure,
     &           maxngrid_r,ngrid_r,grid_r)


        do i=1,ndep_solid
         if(ir_for_stress(1,i).eq.ir_for_stress(2,i).or.
     &      ir_for_stress(2,i).eq.ir_for_stress(3,i)) 
     &     stop 'The upper and lower igrid for stress are the same'
        end do
        do i=1,ndep_fluid
         if(ir_for_pressure(1,i).eq.ir_for_pressure(2,i).or.
     &      ir_for_pressure(2,i).eq.ir_for_pressure(3,i))
     &     stop 'The upper and lower igrid for pressure are the same'
        end do


        ir_CMB=-1;ir_ICB=-1
        do 230 i=1,ngrid_r
           if(dabs(grid_r(i)-r_CMB).lt.1.e-7) ir_CMB=i
           if(dabs(grid_r(i)-r_ICB).lt.1.e-7) ir_ICB=i
  230   continue
        if(ir_CMB.lt.0.or.ir_ICB.lt.0) 
     &            stop 'Not found ir_CMB or ir_ICB'


!################ir_dep_solid and ir_dep_fluid**************
        ir_dep_solid(1:max_ndep)=0
        ir_dep_fluid(1:max_ndep)=0
        r_sta_solid(1:max_ndep)=0.d0
        r_sta_fluid(1:max_ndep)=0.d0

        do 320 i=1,ndep_solid
           depth=depth_solid(i)
           radius=grid_r(ngrid_r)-depth
           r_sta_solid(i)=radius
           ir=1
           do 310 while(grid_r(ir).lt.radius.and.ir<=ngrid_r)
               ir=ir+1
  310      continue
           ir_dep_solid(i)=ir
  320    continue

        do 420 i=1,ndep_fluid
           depth=depth_fluid(i)
           radius=grid_r(ngrid_r)-depth
           r_sta_fluid(i)=radius
           ir=1
           do 410 while(grid_r(ir).lt.radius.and.ir<=ngrid_r)
               ir=ir+1
  410      continue
           ir_dep_fluid(i)=ir
  420    continue

!crossing ICB or CMB
        if(ndep_solid.gt.0.and.ndep_fluid.gt.0) then
           if(grid_r(ngrid_r)-depth_solid(1).ge.grid_r(ir_CMB)) then
              if(dabs(grid_r(ngrid_r)-depth_solid(ndep_solid)
     &                 -grid_r(ir_CMB)).gt.0.01.or.
     &           dabs(grid_r(ngrid_r)-depth_fluid(1)
     &                 -grid_r(ir_CMB)).gt.0.01)
     &               stop 'Error, crossing CMB, but given depth is 
     &                     not acurately located on CMB'
              ir_dep_solid(ndep_solid)=ir_CMB
              ir_dep_fluid(1)=ir_CMB
           else
              if(dabs(grid_r(ngrid_r)-depth_solid(1)
     &                 -grid_r(ir_ICB)).gt.0.01.or.
     &           dabs(grid_r(ngrid_r)-depth_fluid(ndep_fluid)
     &                 -grid_r(ir_ICB)).gt.0.01)
     &               stop 'Error, crossing ICB, but given depth is 
     &                     not acurately located on ICB!!!'
              ir_dep_solid(1)=ir_ICB
              ir_dep_fluid(ndep_fluid)=ir_ICB
           end if
        end if

c        max1_dr=maxval(grid_r(ir_dep_solid(ndep_solid):ir_dep_solid(1)))
c        max2_dr=maxval(grid_r(ir_dep_fluid(ndep_fluid):ir_dep_fluid(1)))
cccccccc   if(max1_dr.gt.ddep.or.max2_dr.gt.ddep) 
ccc     &           stop 'max_dr is larger than ddep'
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine assign_structure
     &	  ( ngrid_r,grid_r,
     &	    n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &	    rho_structure_zone,
     &	    vpv_structure_zone,vph_structure_zone,
     &	    vsv_structure_zone,vsh_structure_zone,
     &	    eta_structure_zone,
     &	    qkappa_structure_zone,qmu_structure_zone,
     &	    grid_rho,grid_Ak,grid_Am,grid_Ck,grid_Cm,
     &	    grid_L,grid_N,grid_Fk,grid_Fm,
     &	    grid_kappa,grid_mu,grid_qkappa,grid_qmu,
     &	    idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &      A_sta,C_sta,L_sta,N_sta,F_sta )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c assigning the input model to the numerical grids
c    required subroutines: error_handling
c    required function: cal_PREM_structure
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer ngrid_r,n_structure_zone
	integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
	real*8 grid_r(*)
	real*8 rmin_structure_zone(*),rmax_structure_zone(*)
	real*8 rho_structure_zone(4,*)
	real*8 vpv_structure_zone(4,*),vph_structure_zone(4,*)
	real*8 vsv_structure_zone(4,*),vsh_structure_zone(4,*)
	real*8 eta_structure_zone(4,*)
	real*8 qkappa_structure_zone(*),qmu_structure_zone(*)
	real*8 grid_rho(2,*)
	real*8 grid_Ak(2,*),grid_Am(2,*),grid_Ck(2,*),grid_Cm(2,*)
	real*8 grid_L(2,*),grid_N(2,*),grid_Fk(2,*),grid_Fm(2,*)
	real*8 grid_kappa(2,*),grid_mu(2,*)
	real*8 grid_qkappa(*),grid_qmu(*)
        real*8 A_sta(*),C_sta(*),L_sta(*),N_sta,F_sta
c other variables
	integer ir,izone,izone1,izone2
	real*8 rho_r(2),vpv_r(2),vph_r(2),vsv_r(2),vsh_r(2),eta_r(2)
	real*8 A_r(2),C_r(2),L_r(2),N_r(2),F_r(2)
	real*8 kappa_r(2),mu_r(2)
	real*8 cal_PREM_structure
c
	do 500 ir=1,ngrid_r-1
	  izone1 = 0
	  izone2 = 0
	  do 150 izone=1,n_structure_zone
	    if ( rmin_structure_zone(izone).le.grid_r(ir) )
     &	      izone1 = izone
  150	  continue
	  do 160 izone=n_structure_zone,1,-1
	    if ( rmax_structure_zone(izone).ge.grid_r(ir+1) )
     &	      izone2 = izone
  160	  continue
	  rho_r(1)
     &	    = cal_PREM_structure( grid_r(ir),
     &	                          rho_structure_zone(1,izone1) )
	  rho_r(2)
     &	    = cal_PREM_structure( grid_r(ir+1),
     &	                          rho_structure_zone(1,izone2) )
	  vpv_r(1)
     &	    = cal_PREM_structure( grid_r(ir),
     &	                          vpv_structure_zone(1,izone1) )
	  vpv_r(2)
     &	    = cal_PREM_structure( grid_r(ir+1),
     &	                          vpv_structure_zone(1,izone2) )
	  vph_r(1)
     &	    = cal_PREM_structure( grid_r(ir),
     &	                          vph_structure_zone(1,izone1) )
	  vph_r(2)
     &	    = cal_PREM_structure( grid_r(ir+1),
     &	                          vph_structure_zone(1,izone2) )
	  vsv_r(1)
     &	    = cal_PREM_structure( grid_r(ir),
     &	                          vsv_structure_zone(1,izone1) )
	  vsv_r(2)
     &	    = cal_PREM_structure( grid_r(ir+1),
     &	                          vsv_structure_zone(1,izone2) )
	  vsh_r(1)
     &	    = cal_PREM_structure( grid_r(ir),
     &	                          vsh_structure_zone(1,izone1) )
	  vsh_r(2)
     &	    = cal_PREM_structure( grid_r(ir+1),
     &	                          vsh_structure_zone(1,izone2) )
	  eta_r(1)
     &	    = cal_PREM_structure( grid_r(ir),
     &	                          eta_structure_zone(1,izone1) )
	  eta_r(2)
     &	    = cal_PREM_structure( grid_r(ir+1),
     &	                          eta_structure_zone(1,izone2) )
	  A_r(1) = rho_r(1) * vph_r(1) * vph_r(1)
	  A_r(2) = rho_r(2) * vph_r(2) * vph_r(2)
	  C_r(1) = rho_r(1) * vpv_r(1) * vpv_r(1)
	  C_r(2) = rho_r(2) * vpv_r(2) * vpv_r(2)
	  L_r(1) = rho_r(1) * vsv_r(1) * vsv_r(1)
	  L_r(2) = rho_r(2) * vsv_r(2) * vsv_r(2)
	  N_r(1) = rho_r(1) * vsh_r(1) * vsh_r(1)
	  N_r(2) = rho_r(2) * vsh_r(2) * vsh_r(2)
	  F_r(1) = eta_r(1) * ( A_r(1) - 2.d0 * L_r(1) )
	  F_r(2) = eta_r(2) * ( A_r(2) - 2.d0 * L_r(2) )

	  kappa_r(1) = ( 4.d0 * A_r(1) + C_r(1)
     &	                   + 4.d0 * F_r(1) - 4.d0 * N_r(1) ) / 9.d0
	  kappa_r(2) = ( 4.d0 * A_r(2) + C_r(2)
     &	                   + 4.d0 * F_r(2) - 4.d0 * N_r(2) ) / 9.d0
	  mu_r(1) = ( A_r(1) + C_r(1) - 2.d0 * F_r(1)
     &	                + 5.d0 * N_r(1) + 6.d0 * L_r(1) ) / 15.d0
	  mu_r(2) = ( A_r(2) + C_r(2) - 2.d0 * F_r(2)
     &	                + 5.d0 * N_r(2) + 6.d0 * L_r(2) ) / 15.d0
	  grid_rho(1,ir) = rho_r(1)
	  grid_rho(2,ir) = rho_r(2)
	  grid_Ak(1,ir) = A_r(1)
     &	                    * kappa_r(1)
     &	                      / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
	  grid_Ak(2,ir) = A_r(2)
     &	                    * kappa_r(2)
     &	                      / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
	  grid_Am(1,ir) = A_r(1)
     &	                    * 4.d0/3.d0*mu_r(1)
     &	                      / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
	  grid_Am(2,ir) = A_r(2)
     &	                    * 4.d0/3.d0*mu_r(2)
     &	                      / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
	  grid_Ck(1,ir) = C_r(1)
     &	                    * kappa_r(1)
     &	                      / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
	  grid_Ck(2,ir) = C_r(2)
     &	                    * kappa_r(2)
     &	                      / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
	  grid_Cm(1,ir) = C_r(1)
     &	                    * 4.d0/3.d0*mu_r(1)
     &	                      / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
	  grid_Cm(2,ir) = C_r(2)
     &	                    * 4.d0/3.d0*mu_r(2)
     &	                      / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
	  grid_L(1,ir) = L_r(1)
	  grid_L(2,ir) = L_r(2)
	  grid_N(1,ir) = N_r(1)
	  grid_N(2,ir) = N_r(2)
	  grid_Fk(1,ir) = F_r(1)
     &	                    * kappa_r(1)
     &	                      / ( kappa_r(1)-2.d0/3.d0*mu_r(1) )
	  grid_Fk(2,ir) = F_r(2)
     &	                    * kappa_r(2)
     &	                      / ( kappa_r(2)-2.d0/3.d0*mu_r(2) )
	  grid_Fm(1,ir) = F_r(1)
     &	                    * (-2.d0/3.d0)*mu_r(1)
     &	                      / ( kappa_r(1)-2.d0/3.d0*mu_r(1) )
	  grid_Fm(2,ir) = F_r(2)
     &	                    * (-2.d0/3.d0)*mu_r(2)
     &	                      / ( kappa_r(2)-2.d0/3.d0*mu_r(2) )
	  grid_kappa(1,ir) = kappa_r(1)
	  grid_kappa(2,ir) = kappa_r(2)
	  grid_mu(1,ir) = mu_r(1)
	  grid_mu(2,ir) = mu_r(2)
	  grid_qkappa(ir)
     &	    = 2.d0/(1.d0/qkappa_structure_zone(izone1)
     &	            +1.d0/qkappa_structure_zone(izone2))
	  grid_qmu(ir)
     &	    = 2.d0/(1.d0/qmu_structure_zone(izone1)
     &	            +1.d0/qmu_structure_zone(izone2))
  500	continue

c --- computing positions of in solution vectors
	idim1_sph = 1
	idim2_sph = ngrid_r
	idim1_tor = 0
	idim2_tor = 0
	do 600 ir=1,ngrid_r-1
	  if ( grid_mu(1,ir)*grid_mu(2,ir).ne.0.d0 ) idim2_tor = ir+1
  600	continue
	if ( idim2_tor.eq.0 ) call error_handling(16)
	do 700 ir=1,idim2_tor-1
	  if ( grid_mu(1,ir)*grid_mu(2,ir).eq.0.d0 ) then
	    idim1_tor = ir
	  endif
  700	continue
	idim1_tor = idim1_tor + 1
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine assign_source
     &	    ( ngrid_r,grid_r,source_r,
     &	      grid_L,grid_Ck,grid_Cm,grid_Fk,grid_Fm,
     &	      grid_qkappa,grid_qmu,grid_mu,
     &	      igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &	      grid_qkappas,grid_qmus,
     &	      idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing intermediate parameters used in assing source parameters
c to the numerical grids
c    required subroutines: error_handling
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer ngrid_r,igrid_rs
	integer idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0
	real*8 grid_r(*),source_r
	real*8 grid_L(2,*),grid_Ck(2,*),grid_Cm(2,*)
	real*8 grid_Fk(2,*),grid_Fm(2,*)
	real*8 grid_qkappa(*),grid_qmu(*),grid_mu(2,*)
	real*8 grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
	real*8 grid_qkappas,grid_qmus
c other variables
	integer i,ipos1,ipos2,ipos3
c
c ---- serching the grid on the source
	igrid_rs = 0
	do 100 i=1,ngrid_r
	  if ( grid_r(i).eq.source_r ) igrid_rs = i
  100	continue
	if ( igrid_rs.eq.0 ) call error_handling(17)
c ---- computing the elastic constants at the source
	if ( igrid_rs.eq.ngrid_r ) then
	  grid_Ls = grid_L(2,igrid_rs-1)
	  grid_Cks = grid_Ck(2,igrid_rs-1)
	  grid_Cms = grid_Cm(2,igrid_rs-1)
	  grid_Fks = grid_Fk(2,igrid_rs-1)
	  grid_Fms = grid_Fm(2,igrid_rs-1)
	  grid_qkappas = grid_qkappa(igrid_rs-1)
	  grid_qmus = grid_qmu(igrid_rs-1)
	else
	  grid_Ls = grid_L(1,igrid_rs)
	  grid_Cks = grid_Ck(1,igrid_rs)
	  grid_Cms = grid_Cm(1,igrid_rs)
	  grid_Fks = grid_Fk(1,igrid_rs)
	  grid_Fms = grid_Fm(1,igrid_rs)
	  grid_qkappas = grid_qkappa(igrid_rs)
	  grid_qmus = grid_qmu(igrid_rs)
	endif
c --- computing positions of non-zero elements in excitation vectors
	ipos1 = 0
	ipos3 = 0
	ipos2 = 0
	if ( igrid_rs.gt.1 ) then
	  do 200 i=1,igrid_rs-1
	    if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
	      ipos1 = ipos1 + 1
	      ipos3 = ipos3 + 1
	      ipos2 = 0
	    else
	      ipos1 = ipos1 + 2
	      ipos3 = ipos3 + 1
	      ipos2 = ipos2 + 1
	    endif
	    if ( ( grid_mu(2,i).eq.0.d0 ).and.
     &	           ( grid_mu(1,i+1).ne.0.d0 )    ) then
	      ipos1 = ipos1 + 1
	      ipos3 = ipos3 + 1
	      ipos2 = 0
	    endif
	    if ( ( grid_mu(2,i).ne.0.d0 ).and.
     &	           ( grid_mu(1,i+1).eq.0.d0 )    ) then
	      ipos1 = ipos1 + 2
	      ipos3 = ipos3 + 1
	      ipos2 = 0
	    endif
  200	  continue
	endif
	idim_rs_sph = ipos1 + 1
	idim_rs_tor = ipos2 + 1
	idim_rs_sph0 = ipos3 + 1
	idim_rs_tor0 = 0
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c WENBO
	subroutine assign_station
     &	    ( ngrid_r,grid_r,grid_mu,
     &        idim_sph_forcut,idim_tor_forcut,
     &        idim_sph_forcut0,idim_tor_forcut0,
     &        minidim_sph_sta,minidim_tor_sta,
     &        minidim_sph_sta0,minidim_tor_sta0,
     &        idim_ir_sph,idim_ir_tor,idim_ir_sph0,
     &        ir_dep_solid,ir_dep_fluid,ndep_solid,
     &        ndep_fluid,idep_forcut,min_ir_dep)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing intermediate parameters used in assing station parameters
c    required subroutines: none
c    required function: error_handling
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
        implicit none
	integer ngrid_r
c	integer idim_station_sph,idim_station_tor,
c     &	        idim_station_sph0,idim_station_tor0
	real*8 grid_r(*),grid_mu(2,*)
c WENBO
        integer idim_sph_forcut,idim_tor_forcut,
     &          idim_sph_forcut0,idim_tor_forcut0,
     &          minidim_sph_sta,minidim_tor_sta,
     &          minidim_sph_sta0,minidim_tor_sta0

        integer idim_ir_sph(*),idim_ir_tor(*),
     &          idim_ir_sph0(*)
        integer ir_dep_solid(*),ir_dep_fluid(*)
        integer ndep_solid,ndep_fluid,idep_forcut
c other variables
	integer ir,idim,i,ipos1,ipos2,ipos3
        integer ir_dep_forcut,min_ir_dep
c
	idim = 0
	do 600 ir=1,ngrid_r-1
	  if ( grid_mu(1,ir)*grid_mu(2,ir).ne.0.d0 ) idim = ir
  600	continue
	if ( idim.eq.0 ) call error_handling(18)
c
c --- computing positions of non-zero elements in excitation vectors
	ipos1 = 0
	ipos3 = 0
	ipos2 = 0
	if ( idim.gt.1 ) then
	  do 200 i=1,idim-1
	    if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
	      ipos1 = ipos1 + 1
	      ipos3 = ipos3 + 1
	      ipos2 = 0
	    else
	      ipos1 = ipos1 + 2
	      ipos3 = ipos3 + 1
	      ipos2 = ipos2 + 1
	    endif
	    if ( ( grid_mu(2,i).eq.0.d0 ).and.
     &	           ( grid_mu(1,i+1).ne.0.d0 )    ) then
	      ipos1 = ipos1 + 1
	      ipos3 = ipos3 + 1
	      ipos2 = 0
	    endif
	    if ( ( grid_mu(2,i).ne.0.d0 ).and.
     &	           ( grid_mu(1,i+1).eq.0.d0 )    ) then
	      ipos1 = ipos1 + 2
	      ipos3 = ipos3 + 1
	      ipos2 = 0
	    endif
            idim_ir_sph(i)=ipos1-1
            idim_ir_tor(i)=ipos2
            idim_ir_sph0(i)=ipos3

  200	  continue
	endif
c	idim_station_sph = ipos1 + 1 + 2
c	idim_station_tor = ipos2 + 1 + 1
c	idim_station_sph0 = ipos3 + 1 + 1
c	idim_station_tor0 = 0
        idim_ir_sph(ngrid_r-1)=ipos1+1
        idim_ir_tor(ngrid_r-1)=ipos2+1
        idim_ir_sph0(ngrid_r-1)=ipos3+1
        idim_ir_sph(ngrid_r)=ipos1+1+2
        idim_ir_tor(ngrid_r)=ipos2+1+1
        idim_ir_sph0(ngrid_r)=ipos3+1+1
        if(ndep_fluid.eq.0) then
          ir_dep_forcut=ir_dep_solid(idep_forcut)
          min_ir_dep=ir_dep_solid(ndep_solid)
        else if(ndep_solid.eq.0) then
          ir_dep_forcut=ir_dep_fluid(idep_forcut)
          min_ir_dep=ir_dep_fluid(ndep_fluid)
        else
          if(ir_dep_fluid(1).lt.ir_dep_solid(1)) then
            min_ir_dep=ir_dep_fluid(ndep_fluid)
          else
            min_ir_dep=ir_dep_solid(ndep_solid)
          end if
          if(idep_forcut.le.ndep_solid) then
              ir_dep_forcut=ir_dep_solid(idep_forcut)
          else
              ir_dep_forcut=ir_dep_fluid(idep_forcut-ndep_solid)
          end if

        end  if

        idim_sph_forcut=idim_ir_sph(ir_dep_forcut)
        idim_tor_forcut=idim_ir_tor(ir_dep_forcut)
        idim_sph_forcut0=idim_ir_sph0(ir_dep_forcut)
        
        minidim_sph_sta=idim_ir_sph(min_ir_dep)
        minidim_tor_sta=idim_ir_tor(min_ir_dep)
        minidim_sph_sta0=idim_ir_sph0(min_ir_dep)

c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 function cal_PREM_structure( r,param )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 r,param(4)
	real*8 rmax,a
c
	rmax = 6371.d0
	a = r / rmax
	cal_PREM_structure
     &	  = param(1)
     &	    + param(2) * a
     &	    + param(3) * a * a
     &	    + param(4) * a * a * a
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine write_spectrum_file
     &	    ( i_frequency,n_station,spectrum_file,
     &	      station_displacement )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c writing the displacement at each station to the spectrum files.
c   required subroutine: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer i_frequency,n_station
	character*80 spectrum_file(*)
	complex*16 station_displacement(3,*)
c other variables
	integer i_station
c
	if ( i_frequency.eq.0 ) then
	  do 100 i_station=1,n_station
	    open( unit=11,file=spectrum_file(i_station),
     &	          status='unknown' )
	    write(11,*) i_frequency,
     &	                station_displacement(1,i_station),
     &	                station_displacement(2,i_station),
     &	                station_displacement(3,i_station)
	    close(11)
  100	  continue
	else
	  do 110 i_station=1,n_station
	    open( unit=11,file=spectrum_file(i_station),
     &	          access='append',status='old' )
	    write(11,*) i_frequency,
     &	                station_displacement(1,i_station),
     &	                station_displacement(2,i_station),
     &	                station_displacement(3,i_station)
	    close(11)
  110	  continue
	endif
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine distaz(the,phe,ths,phs,ns,dist,az,baz,xdeg,nerr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*=====================================================================
* PURPOSE:  To compute the distance and azimuth between locations.
*=====================================================================
* INPUT ARGUMENTS:
*    THE:     Event latitude in decimal degrees, North positive. [r]
*    PHE:     Event longitude, East positive. [r]
*    THS:     Array of station latitudes. [r]
*    PHS:     Array of station longitudes. [r]
*    NS:      Length of THS and PHS. [i]
*=====================================================================
* OUTPUT ARGUMENTS:
*    DIST:    Array of epicentral distances in km. [r]
*    AZ:      Array of azimuths in degrees. [r]
*    BAZ:     Array of back azimuths. [r]
*    XDEG:    Array of great circle arc lengths. [r]
*    NERR:    Error flag:
*             =    0   No error.
*             = 0904   Calculation failed internal consistency checks.
*=====================================================================
* MODULE/LEVEL:  DFM/4
*=====================================================================
* GLOBAL INPUT:
*    MACH:
*=====================================================================
* SUBROUTINES CALLED:
*    SACLIB:  SETMSG, APCMSG
*=====================================================================
* LOCAL VARIABLES:
*=====================================================================
* KNOWN ERRORS:
* - Problem with equation for distance. See discussion below.
*=====================================================================

	real pi,todeg,torad
	parameter (PI=3.141592654)
	parameter (TODEG=57.29577950)
	parameter (TORAD=1./TODEG)
c
	integer ns,nerr,i
	real ths(ns), phs(ns)
	real dist(ns), az(ns), baz(ns), xdeg(ns)
	logical laz,lbaz,ldist,lxdeg
	real the,phe,ec2,onemec2,eps,temp,therad,pherad,thg
	real d,e,f,c,a,b,g,h,thsrad,phsrad
	real d1,e1,f1,c1,a1,b1,g1,h1,sc,sd,ss,t1,p1,t2,p2,el
	real costhi,costhk,sinthi,sinthk,tanthi,tanthk,al,dl
	real a12top,a12bot,a12,cosa12,sina12,e2,e3,c0,c2,c4,v1,v2,z1,z2,x2,y2
	real e1p1,sqrte1p1,u1bot,u1,u2top,u2bot,u2,b0,du,pdist
	real rad,fl,twopideg,degtokm
	real c00,c01,c02,c03,c21,c22,c23,c42,c43
* PROCEDURE:

* - Calculations are based upon the reference spheroid of 1968 and
*   are defined by the major radius (RAD) and the flattening (FL).

      data rad/6378.160/,fl/0.00335293/
      data twopideg/360./
      data c00,c01,c02,c03/1.,0.25,-4.6875e-02,1.953125e-02/
      data c21,c22,c23/-0.125,3.125e-02,-1.46484375e-02/
      data c42,c43/-3.90625e-03,2.9296875e-03/
      data degtokm/111.3199/

* - Initialize.

      nerr=0
      ec2=2.*fl-fl*fl
      onemec2=1.-ec2
      eps=1.+ec2/onemec2

* - Check which output items are required.

      laz=.true.
c      if(az(1).lt.0.)laz=.false.
      lbaz=.true.
c      if(baz(1).lt.0.)lbaz=.false.
      ldist=.true.
c      if(dist(1).lt.0.)ldist=.false.
      lxdeg=.true.
c      if(xdeg(1).lt.0.)lxdeg=.false.

* - Convert event location to radians.
*   (Equations are unstable for latidudes of exactly 0 degrees.)

      temp=the
      if(temp.eq.0.)temp=1.0e-08
      therad=torad*temp
      pherad=torad*phe

* - Must convert from geographic to geocentric coordinates in order
*   to use the spherical trig equations.  This requires a latitude
*   correction given by: 1-EC2=1-2*FL+FL*FL

      thg=atan(onemec2*tan(therad))
      d=sin(pherad)
      e=-cos(pherad)
      f=-cos(thg)
      c=sin(thg)
      a= f*e
      b=-f*d
      g=-c*e
      h=c*d

* - Loop on stations:

      do 5000 i=1,ns

* -- Convert to radians.
        temp=ths(i)
        if(temp.eq.0.)temp=1.0e-08
        thsrad=torad*temp
        phsrad=torad*phs(i)

* -- Calculate some trig constants.
        thg=atan(onemec2*tan(thsrad))
        d1=sin(phsrad)
        e1=-cos(phsrad)
        f1=-cos(thg)
        c1=sin(thg)
        a1=f1*e1
        b1=-f1*d1
        g1=-c1*e1
        h1=c1*d1
        sc=a*a1+b*b1+c*c1

* - Spherical trig relationships used to compute angles.

        if(lxdeg)then
          sd=0.5*sqrt(((a-a1)**2+(b-b1)**2+(c-c1)**2)*((a+a1)**2
     #       +(b+b1)**2+(c+c1)**2))
          xdeg(i)=atan2(sd,sc)*todeg
          if(xdeg(i).lt.0.)xdeg(i)=xdeg(i)+twopideg
        endif
        if(laz)then
          ss = ((a1-d)**2+(b1-e)**2+c1**2-2.)
          sc = ((a1-g)**2+(b1-h)**2+(c1-f)**2-2.)
          az(i)=atan2(ss,sc)*todeg
          if(az(i).lt.0.)az(i)=az(i)+twopideg
        endif
        if(lbaz)then
          ss=((a-d1)**2+(b-e1)**2+c**2-2.)
          sc=((a-g1)**2+(b-h1)**2+(c-f1)**2-2.)
          baz(i)=atan2(ss,sc)*todeg
          if(baz(i).lt.0.)baz(i)=baz(i)+twopideg
        endif

* - Now compute the distance between the two points using Rudoe's
*   formula given in GEODESY, section 2.15(b).
*   (There is some numerical problem with the following formulae.
*   If the station is in the southern hemisphere and the event in
*   in the northern, these equations give the longer, not the
*   shorter distance between the two locations.  Since the equations
*   are fairly messy, the simplist solution is to reverse the
*   meanings of the two locations for this case.)
        if(ldist)then
          if(thsrad.gt.0.)then
            t1=thsrad
            p1=phsrad
            t2=therad
            p2=pherad
          else
            t1=therad
            p1=pherad
            t2=thsrad
            p2=phsrad
          endif
          el=ec2/onemec2
          e1=1.+el
          costhi=cos(t1)
          costhk=cos(t2)
          sinthi=sin(t1)
          sinthk=sin(t2)
          tanthi=sinthi/costhi
          tanthk=sinthk/costhk
          al=tanthi/(e1*tanthk)+
     #       ec2*sqrt((e1+tanthi**2)/(e1+tanthk**2))
          dl=p1-p2
          a12top=sin(dl)
          a12bot=(al-cos(dl))*sinthk
          a12=atan2(a12top,a12bot)
          cosa12=cos(a12)
          sina12=sin(a12)
          e1=el*((costhk*cosa12)**2+sinthk**2)
          e2=e1*e1
          e3=e1*e2
          c0=c00+c01*e1+c02*e2+c03*e3
          c2=c21*e1+c22*e2+c23*e3
          c4=c42*e2+c43*e3
          v1=rad/sqrt(1.-ec2*sinthk**2)
          v2=rad/sqrt(1.-ec2*sinthi**2)
          z1=v1*(1.-ec2)*sinthk
          z2=v2*(1.-ec2)*sinthi
          x2=v2*costhi*cos(dl)
          y2=v2*costhi*sin(dl)
          e1p1=e1+1.
          sqrte1p1=sqrt(e1p1)
          u1bot=sqrte1p1*cosa12
          u1=atan2(tanthk,u1bot)
          u2top=v1*sinthk+e1p1*(z2-z1)
          u2bot=sqrte1p1*(x2*cosa12-y2*sinthk*sina12)
          u2=atan2(u2top,u2bot)
          b0=v1*sqrt(1.+el*(costhk*cosa12)**2)/e1p1
          du=u2 -u1
          pdist=b0*(c2*(sin(2.*u2)-sin(2.*u1))+
     #       c4*(sin(4.*u2)-sin(4.*u1)))
          dist(i)=abs(b0*c0*du+pdist)
          if(lxdeg .and. (abs(dist(i)-degtokm*xdeg(i))).gt.100.)then
            nerr=0904
c            call setmsg('ERROR',nerr)
c            call apimsg(i)
          endif
        endif
 5000   continue

 8888 return

*=====================================================================
* MODIFICATION HISTORY:
*    830603:  Fixed bug with negative station latiudes.
*    810000:  Original version.
*=====================================================================

      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine error_handling(id)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c error handling.
c    required subroutines: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer id
	integer lnblnk
	character*72 message(99)
c
	message(1) = 'ngrid_r is too large (parameter_input).'
	message(2) = 'n_structure_zone is too large (parameter_input).'
	message(3) = 'lmax is too large (parameter_input).'
	message(5) = 'n_station is too large (parameter_input).'
	message(11) = 'ngrid_r is too large
     & (grid_generation).'
	message(16) = 'Something is wrong (assign_structure).'
	message(17) = 'Something is wrong (assign_source).'
	message(18) = 'Something is wrong (assign_station).'
	message(40) = 'Bad arguments (comp_vecsph_station)'
	message(41) = 'Bad arguments (comp_vecsph_all)'
	message(42) = 'Bad arguments (plgndr)'
	message(51) = 'Bad arguments (comp_excitation)'
	message(51) = 'Bad arguments (comp_excitation0)'
	message(53) = 'Bad arguments (comp_wavefield)'
	message(54) = 'Bad arguments (comp_wavefield0)'
	message(55) = 'Bad arguments (comp_displacement_station)'
	message(56) = 'Bad arguments (comp_displacement_station0)'
c
	write(6,*) '******************** ERROR ********************'
	write(6,*) message(id)(1:lnblnk(message(id)))
	write(6,*) '***********************************************'
	stop
c
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine check_amp_significance
     &	  ( maxngrid_r,l,nl_check_amp,amp_max,i_significance,idim0,
     &	    init_npos_sph,init_npos_tor,ifreq,
     &	    idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &	    grid_mu,idim_station_sph,idim_station_tor,
     &	    whole_vector_sph,whole_vector_tor,work_vector,
     &      min_ir_dep)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
        implicit none
        integer ifreq
	integer maxngrid_r,l,nl_check_amp,i_significance,idim0
	integer init_npos_sph,init_npos_tor
	integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
	integer idim_station_sph,idim_station_tor
	real*8 amp_max,grid_mu(2,*)
	complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
	complex*16 whole_vector_tor(maxngrid_r,-2:2)
	real*8 work_vector(maxngrid_r,2)

c       WENBO
        integer min_ir_dep
c other variables
	integer i,ipos1,ipos2,itmp,idim0_sph,idim0_tor
	real*8 amp,eps,amaxr_sph,amaxr_tor,fac,tinyvalue
c
	data eps / 1.d-12 /
c The precision of real*8 type is 2.23d-308, thus we set tiny value as
c   2.23d-308/1.d-12*10=2.23d-295
        data tinyvalue / 1.d-295 /
c
	if ( mod(l,nl_check_amp).eq.0 ) then
	  if ( l.eq.0 ) then
	    amp_max = cdabs( whole_vector_sph(idim_station_sph,0) )**2
            if(amp_max.lt.tinyvalue) amp_max=tinyvalue
	    i_significance = 1
c
	    ipos1 = 0
	    amaxr_sph = 0.d0
	    amaxr_tor = 0.d0
	    call init_complex_array( maxngrid_r,work_vector )
	    if ( idim2_sph.gt.idim1_sph ) then
c introducing 'fac' to avoid underflow exceptions
	      fac = 1.d100
	      do 100 i=idim1_sph,idim2_sph-1
	        if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
	          ipos1 = ipos1 + 1
	        else
	          ipos1 = ipos1 + 1
	          work_vector(i,1)
     &	            = cdabs( whole_vector_sph(ipos1,0)*fac )**2
	          if ( work_vector(i,1).gt.amaxr_sph )
     &	            amaxr_sph = work_vector(i,1)
	        endif
	        if ( ( grid_mu(2,i).eq.0.d0 ).and.
     &	               ( grid_mu(1,i+1).ne.0.d0 )    ) then
	          ipos1 = ipos1 + 1
	        endif
	        if ( ( grid_mu(2,i).ne.0.d0 ).and.
     &	               ( grid_mu(1,i+1).eq.0.d0 )    ) then
	          ipos1 = ipos1 + 1
	        endif
  100	      continue
	      itmp = idim2_sph
	      do 110 i=idim2_sph-1,idim1_sph,-1
	        if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
	          if ( itmp.eq.i+1 ) itmp = i
	        else
c	          if ( work_vector(i,1).gt.eps*amaxr_sph ) itmp=i
	          if ( work_vector(i,1).ge.eps*amaxr_sph ) itmp=i
	        endif
  110	      continue
	      if ( itmp.ne.idim2_sph ) idim0_sph = itmp
	      if ( idim0_sph.gt.2 ) idim0 = idim0_sph
	    endif
c
	  else
	    amp =
     &	        cdabs( whole_vector_sph(idim_station_sph,  -2) )**2
     &	      + cdabs( whole_vector_sph(idim_station_sph+1,-2) )**2
     &	      + cdabs( whole_vector_tor(idim_station_tor,  -2) )**2
     &	      + cdabs( whole_vector_sph(idim_station_sph,  -1) )**2
     &	      + cdabs( whole_vector_sph(idim_station_sph+1,-1) )**2
     &	      + cdabs( whole_vector_tor(idim_station_tor,  -1) )**2
     &	      + cdabs( whole_vector_sph(idim_station_sph,   0) )**2
     &	      + cdabs( whole_vector_sph(idim_station_sph+1, 0) )**2
     &	      + cdabs( whole_vector_sph(idim_station_sph,   1) )**2
     &	      + cdabs( whole_vector_sph(idim_station_sph+1, 1) )**2
     &	      + cdabs( whole_vector_tor(idim_station_tor,   1) )**2
     &	      + cdabs( whole_vector_sph(idim_station_sph,   2) )**2
     &	      + cdabs( whole_vector_sph(idim_station_sph+1, 2) )**2
     &	      + cdabs( whole_vector_tor(idim_station_tor,   2) )**2
	    if ( amp.gt.amp_max ) amp_max = amp

cWENBO ------   0: do one more round. Take the averaged solution
c                  (between l_sig0 and l_sig0+nl_check_amp) to 
c                  reduce the osilation around source depth.
cWENBO ------  -1: calculation completed
            if(i_significance.eq.0) then
c               write(6,*) 'i_significance=-1'
               i_significance=-1
	    else if ( amp.lt.eps*amp_max ) then
               i_significance = 0
            end if

cWENBO*************************************************************



c
	    ipos1 = 0
	    ipos2 = 0
	    amaxr_sph = 0.d0
	    amaxr_tor = 0.d0
	    call init_complex_array( maxngrid_r,work_vector )
	    ipos1 = 0
	    ipos2 = 0
	    if ( idim2_sph.gt.idim1_sph ) then
c introducing 'fac' to avoid underflow exceptions
	      fac = 1.d100
	      do 200 i=idim1_sph,idim2_sph-1
	        if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
	          ipos1 = ipos1 + 1
	        else
	          ipos1 = ipos1 + 2
	          work_vector(i,1)
     &	            =   cdabs( whole_vector_sph(ipos1,  -2)*fac )**2
     &	              + cdabs( whole_vector_sph(ipos1+1,-2)*fac )**2
     &	              + cdabs( whole_vector_sph(ipos1,  -1)*fac )**2
     &	              + cdabs( whole_vector_sph(ipos1+1,-1)*fac )**2
     &	              + cdabs( whole_vector_sph(ipos1,   0)*fac )**2
     &	              + cdabs( whole_vector_sph(ipos1+1, 0)*fac )**2
     &	              + cdabs( whole_vector_sph(ipos1,   1)*fac )**2
     &	              + cdabs( whole_vector_sph(ipos1+1, 1)*fac )**2
     &	              + cdabs( whole_vector_sph(ipos1,   2)*fac )**2
     &	              + cdabs( whole_vector_sph(ipos1+1, 2)*fac )**2
	          if ( work_vector(i,1).gt.amaxr_sph )
     &	            amaxr_sph = work_vector(i,1)
	        endif
	        if ( ( grid_mu(2,i).eq.0.d0 ).and.
     &	               ( grid_mu(1,i+1).ne.0.d0 )    ) then
	          ipos1 = ipos1 + 1
	        endif
	        if ( ( grid_mu(2,i).ne.0.d0 ).and.
     &	               ( grid_mu(1,i+1).eq.0.d0 )    ) then
	          ipos1 = ipos1 + 2
	        endif
  200	      continue
c
	      do 210 i=idim1_tor,idim2_tor-1
	        ipos2 = ipos2 + 1
	        work_vector(i,2)
     &	          =   cdabs( whole_vector_tor(ipos2,-2)*fac )**2
     &	            + cdabs( whole_vector_tor(ipos2,-1)*fac )**2
     &	            + cdabs( whole_vector_tor(ipos2, 1)*fac )**2
     &	            + cdabs( whole_vector_tor(ipos2, 2)*fac )**2
	        if ( work_vector(i,2).gt.amaxr_tor )
     &	          amaxr_tor = work_vector(i,2)
  210	      continue
c
	      itmp = idim2_sph
	      do 220 i=idim2_sph-1,idim1_sph,-1
	        if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
	          if ( itmp.eq.i+1 ) itmp = i
	        else
c	          if ( work_vector(i,1).gt.eps*amaxr_sph ) itmp=i
	          if ( work_vector(i,1).ge.eps*amaxr_sph ) itmp=i
	        endif
  220	      continue

ccccc WENBO  Even if the coefficient amplitude at station grids are
cccc         under threshold, we still solve it.
              if(itmp.gt.min_ir_dep) itmp=min_ir_dep

ccccc WENBO----renew idim0 to take more nodes into calculations
	      if ( itmp.ne.idim2_sph ) then
	        idim0_sph = itmp
	      else
	        idim0_sph = idim1_sph
	      endif
	      itmp = idim2_tor
	      do 230 i=idim2_tor-1,idim1_tor,-1
c	        if ( work_vector(i,2).gt.eps*amaxr_tor ) itmp=i
	        if ( work_vector(i,2).ge.eps*amaxr_tor ) itmp=i
  230	      continue

cccc WEENBO Even if the coefficient amplitude at station grids are
cccc         under threshold, we still solve it.
              if(itmp.gt.min_ir_dep) itmp=min_ir_dep

	      if ( itmp.ne.idim2_tor ) then
	        idim0_tor = itmp
	      else
	        idim0_tor = idim1_tor
	      endif
	      if ( amaxr_tor.ne.0.d0 ) then
	        if ( ( idim0.lt.min0(idim0_sph,idim0_tor) ).and.
     &	             ( min0(idim0_sph,idim0_tor).gt.2 ) )
     &	          idim0 = min0(idim0_sph,idim0_tor)
	      else
	        if ( idim0_sph.gt.2 )
     &	          idim0 = idim0_sph
	      endif
	    endif
	  endif
c
	  ipos1 = 0
	  ipos2 = 0
	  if ( idim0.gt.1 ) then
	    do 300 i=1,idim0-1
	      if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
	        ipos1 = ipos1 + 1
	        ipos2 = 0
	      else
	        ipos1 = ipos1 + 2
	        ipos2 = ipos2 + 1
	      endif
	      if ( ( grid_mu(2,i).eq.0.d0 ).and.
     &	             ( grid_mu(1,i+1).ne.0.d0 )    ) then
	        ipos1 = ipos1 + 1
	        ipos2 = 0
	      endif
	      if ( ( grid_mu(2,i).ne.0.d0 ).and.
     &	             ( grid_mu(1,i+1).eq.0.d0 )    ) then
	        ipos1 = ipos1 + 2
	        ipos2 = 0
	      endif
  300	    continue
	  endif
	  init_npos_sph = ipos1 + 1
	  if ( idim0.ge.idim1_tor ) then
	    init_npos_tor = ipos2 + 1
	  else
	    init_npos_tor = 1
	  endif
c
	endif
c
	end
