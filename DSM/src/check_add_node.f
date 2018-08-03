        subroutine addnodes_besides_topbot(top_fluid,bot_fluid,
     &           nstructure_zone,rmax_structure_zone,
     &           rmin_structure_zone,izone_idep_solid,
     &           izone_idep_fluid,depth_tolerence,
     &           depth_solid,depth_fluid,ndep_solid,
     &           ndep_fluid,depth_for_stress,depth_for_pressure,
     &           ir_for_stress,ir_for_pressure,r_top_stress,
     &           r_bot_stress,r_top_pressure,r_bot_pressure,
     &           ir_top_stress,ir_bot_stress,
     &           ir_top_pressure,ir_bot_pressure,
     &           maxngrid_r,ngrid_r,grid_r)
        implicit none
        integer ndep_solid,ndep_fluid,nstructure_zone
        integer ngrid_r,maxngrid_r
        integer izone_idep_solid(*),izone_idep_fluid(*)
        integer ir_for_stress(3,*),ir_for_pressure(3,*)
        integer ir_top_stress(3),ir_bot_stress(3)
        integer ir_top_pressure(3),ir_bot_pressure(3)
        real*8 depth_tolerence
        real*8 depth_for_stress(3,*),depth_for_pressure(3,*)
        real*8 depth_solid(*),depth_fluid(*)
        real*8 r_top_stress(3),r_bot_stress(3)
        real*8 r_top_pressure(3),r_bot_pressure(3)
        real*8 rmax_structure_zone(*),rmin_structure_zone(*)
        real*8 grid_r(*)
        integer top_fluid,bot_fluid
!local parameters
        integer idep,i,i_updown,ir_add_node
        real*8 r_freesurf,r_add_node,r_top,r_bot
        real*8 rmin_thiszone,rmax_thiszone
        
        r_freesurf=rmax_structure_zone(nstructure_zone)
        if(ndep_solid.gt.0) then
         do idep=ndep_solid,1,-1
          rmax_thiszone=rmax_structure_zone(izone_idep_solid(idep))
          rmin_thiszone=rmin_structure_zone(izone_idep_solid(idep))
          do i_updown=1,3 !1 for down, 2 for mid, and 3 for up

           r_add_node=r_freesurf-depth_for_stress(i_updown,idep)
           call check_add_node(r_add_node,ir_add_node,grid_r,
     &              ngrid_r,maxngrid_r,depth_tolerence,
     &              ndep_solid,ndep_fluid,
     &              ir_for_stress,ir_for_pressure,ir_top_stress,
     &              ir_bot_stress,ir_top_pressure,ir_bot_pressure)

           if(grid_r(ir_add_node).gt.rmax_thiszone.or.
     &        grid_r(ir_add_node).lt.rmin_thiszone) then
c                 print *,'ir_dep_stress_find',idep,grid_r(ir_add_node),
c     &              depth_for_stress(i_updown,idep)
                 stop 'Error, r_for_stress is not within the 
     &                  izone where the input depth(solid) is'
           end if
           ir_for_stress(i_updown,idep)=ir_add_node
          end do !do i_updown
         end do ! do idep
        end if

        if(ndep_fluid.gt.0) then
         do idep=ndep_fluid,1,-1
          rmax_thiszone=rmax_structure_zone(izone_idep_fluid(idep))
          rmin_thiszone=rmin_structure_zone(izone_idep_fluid(idep))
          do i_updown=1,3
           r_add_node=r_freesurf-depth_for_pressure(i_updown,idep)
           call check_add_node(r_add_node,ir_add_node,grid_r,
     &             ngrid_r,maxngrid_r,depth_tolerence,
     &             ndep_solid,ndep_fluid,
     &             ir_for_stress,ir_for_pressure,ir_top_stress,
     &             ir_bot_stress,ir_top_pressure,ir_bot_pressure)
           if(grid_r(ir_add_node).gt.rmax_thiszone.or.
     &        grid_r(ir_add_node).lt.rmin_thiszone) then
c                 print *,'ir_dep_stress_find',idep,grid_r(ir_add_node)
                 stop 'Error, r_for_stress is not within the 
     &                  izone where the input depth(fluid) is'
           end if
           ir_for_pressure(i_updown,idep)=ir_add_node
          end do !do i_updown
         end do ! do idep
        end if

        
        if(top_fluid.eq.0) then
          r_top=r_freesurf-depth_solid(1)
          rmax_thiszone=rmax_structure_zone(izone_idep_solid(1))
          if(rmax_thiszone-r_top.lt.2.d0*depth_tolerence) then 
            r_top_stress(1)=r_top
            r_top_stress(2)=r_top-1.5*depth_tolerence
            r_top_stress(3)=r_top-3.0*depth_tolerence
          else
            r_top_stress(1)=r_top+1.5*depth_tolerence
            r_top_stress(2)=r_top
            r_top_stress(3)=r_top-1.5*depth_tolerence
          end if

          do i=1,3
            call check_add_node(r_top_stress(i),ir_top_stress(i),
     &             grid_r,ngrid_r,maxngrid_r,depth_tolerence,
     &             ndep_solid,ndep_fluid,
     &             ir_for_stress,ir_for_pressure,ir_top_stress,
     &             ir_bot_stress,ir_top_pressure,ir_bot_pressure)
          end do
        else
          r_top=r_freesurf-depth_fluid(1)
          rmax_thiszone=rmax_structure_zone(izone_idep_fluid(1))
          if(rmax_thiszone-r_top.lt.2.d0*depth_tolerence) then
            r_top_pressure(1)=r_top
            r_top_pressure(2)=r_top-1.5*depth_tolerence
            r_top_pressure(3)=r_top-3.0*depth_tolerence
          else
            r_top_pressure(1)=r_top+1.5*depth_tolerence
            r_top_pressure(2)=r_top
            r_top_pressure(3)=r_top-1.5*depth_tolerence
          end if

          do i=1,3
          call check_add_node(r_top_pressure(i),ir_top_pressure(i),
     &             grid_r,ngrid_r,maxngrid_r,depth_tolerence,
     &             ndep_solid,ndep_fluid,
     &             ir_for_stress,ir_for_pressure,ir_top_stress,
     &             ir_bot_stress,ir_top_pressure,ir_bot_pressure)
          end do
        end if

        if(bot_fluid.eq.0) then
          r_bot=r_freesurf-depth_solid(ndep_solid)
          rmin_thiszone=
     &           rmin_structure_zone(izone_idep_solid(ndep_solid))
          if(r_bot-rmin_thiszone.lt.2.d0*depth_tolerence) then
            r_bot_stress(1)=r_bot+3.0*depth_tolerence
            r_bot_stress(2)=r_bot+1.5*depth_tolerence
            r_bot_stress(3)=r_bot
          else
            r_bot_stress(1)=r_bot+1.5*depth_tolerence
            r_bot_stress(2)=r_bot
            r_bot_stress(3)=r_bot-1.5*depth_tolerence
          end if

          do i=1,3
          call check_add_node(r_bot_stress(i),ir_bot_stress(i),
     &             grid_r,ngrid_r,maxngrid_r,depth_tolerence,
     &             ndep_solid,ndep_fluid,
     &             ir_for_stress,ir_for_pressure,ir_top_stress,
     &             ir_bot_stress,ir_top_pressure,ir_bot_pressure)
          end do
        else
          r_bot=r_freesurf-depth_fluid(ndep_fluid)
          rmin_thiszone=
     &      rmin_structure_zone(izone_idep_fluid(ndep_fluid))
          if(r_bot-rmin_thiszone.lt.2.d0*depth_tolerence) then
            r_bot_pressure(1)=r_bot
            r_bot_pressure(2)=r_bot-1.5*depth_tolerence
            r_bot_pressure(3)=r_bot-3.0*depth_tolerence
          else
            r_bot_pressure(1)=r_bot+1.5*depth_tolerence
            r_bot_pressure(2)=r_bot
            r_bot_pressure(3)=r_bot-1.5*depth_tolerence
          end if

          do i=1,3
          call check_add_node(r_bot_pressure(i),ir_bot_pressure(i),
     &             grid_r,ngrid_r,maxngrid_r,depth_tolerence,
     &             ndep_solid,ndep_fluid,
     &             ir_for_stress,ir_for_pressure,ir_top_stress,
     &             ir_bot_stress,ir_top_pressure,ir_bot_pressure)
          end do
        end if


        end subroutine


        subroutine check_add_node(r_input,ir_return,grid_r,
     &             ngrid_r,maxngrid_r,depth_tolerence,
     &             ndep_solid,ndep_fluid,
     &             ir_for_stress,ir_for_pressure,ir_top_stress,
     &             ir_bot_stress,ir_top_pressure,ir_bot_pressure)
        implicit none
        integer ir_return,ngrid_r,ndep_solid,ndep_fluid,maxngrid_r
        integer ir_for_stress(3,*),ir_for_pressure(3,*)
        integer ir_top_stress(3),ir_bot_stress(3)
        integer ir_top_pressure(3),ir_bot_pressure(3)
        real*8 r_input,depth_tolerence
        real*8 grid_r(*)

!local parameters
        integer i,ir_find,ir_temp,idep
        logical add_node
!find ir_find
        do i=1,ngrid_r
          if(r_input.ge.grid_r(i)) ir_find=i
        end do

        if(ir_find.eq.ngrid_r.and.
     &      (r_input-grid_r(ngrid_r)).gt.depth_tolerence) then
c           print *,'r_input',r_input,grid_r(ngrid_r)
           stop 'r_input is too large' 
        else if (dabs(r_input-grid_r(ir_find)).le.depth_tolerence.and.
     &          ir_find.ge.2) then
           add_node=.false.
        else if (ir_find.lt.ngrid_r.and.
     &         dabs(r_input-grid_r(ir_find+1)).le.depth_tolerence.and.
     &          ir_find.ge.2) then
           add_node=.false.
           ir_find=ir_find+1
        else if(ir_find.gt.2.and.ir_find.lt.ngrid_r) then
           add_node=.true.
           ir_find=ir_find+1
        else
           stop 'Can not find ir_find!'  
        end if



        if(add_node.and.ngrid_r+1.gt.maxngrid_r) 
     &       stop 'Error,ngrid_r>maxngrid_r' 
!move grid_r and change corresponding index
        if(add_node) then
          do ir_temp=ngrid_r,ir_find,-1
            grid_r(ir_temp+1)=grid_r(ir_temp)
          end do
          grid_r(ir_find)=r_input
          ngrid_r=ngrid_r+1
          do idep=1,ndep_solid
            if(ir_for_stress(1,idep).ge.ir_find) then
               ir_for_stress(1,idep)=ir_for_stress(1,idep)+1
            end if
            if(ir_for_stress(2,idep).ge.ir_find) then
               ir_for_stress(2,idep)=ir_for_stress(2,idep)+1
            end if
            if(ir_for_stress(3,idep).ge.ir_find) then
               ir_for_stress(3,idep)=ir_for_stress(3,idep)+1
            end if
          end do
          do idep=1,ndep_fluid
            if(ir_for_pressure(1,idep).ge.ir_find) then
               ir_for_pressure(1,idep)=ir_for_pressure(1,idep)+1
            end if
            if(ir_for_pressure(2,idep).ge.ir_find) then
               ir_for_pressure(2,idep)=ir_for_pressure(2,idep)+1
            end if
            if(ir_for_pressure(3,idep).ge.ir_find) then
               ir_for_pressure(3,idep)=ir_for_pressure(3,idep)+1
            end if
          end do

          do i=1,3
           if(ir_top_stress(i).ge.ir_find) 
     &         ir_top_stress(i)=ir_top_stress(i)+1
           if(ir_bot_stress(i).ge.ir_find) 
     &         ir_bot_stress(i)=ir_bot_stress(i)+1
           if(ir_top_pressure(i).ge.ir_find) 
     &         ir_top_pressure(i)=ir_top_pressure(i)+1
           if(ir_bot_pressure(i).ge.ir_find) 
     &         ir_bot_pressure(i)=ir_bot_pressure(i)+1
          end do
        end if

        ir_return=ir_find

        end subroutine
