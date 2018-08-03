        subroutine write_disp_velo_stress(i_frequency,
     &     time_series_length,omega_imag,
     &     ngrid_r,n_station_solid,n_station_fluid,
     &     r_sta_solid,r_sta_fluid,
     &     ir_dep_solid,max_ndep,maxn_structure_zone,
     &     izone_idep_solid,izone_idep_fluid,
     &     n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &     rho_structure_zone,vpv_structure_zone,vph_structure_zone,
     &     vsv_structure_zone,vsh_structure_zone,eta_structure_zone,
     &     save_velo,ndep_solid,ntheta,
     &     max_nstation,station_disp_solid,
     &     station_disp_fluid,pressure,
     &     epsilon11,epsilon22,epsilon33,
     &     epsilon12,epsilon13,epsilon23,
     &     source_type,top_fluid,bot_fluid,lmax_computed,
     &     lmax,coef_c,coef_dcdr,Atop,Ctop,Ftop,Ltop,Ntop,
     &     Abot,Cbot,Fbot,Lbot,Nbot)
c   input/output
        implicit none
        integer i_frequency,ngrid_r,n_station_solid,
     &          n_station_fluid,max_nstation,ndep_solid,
     &          ntheta,max_ndep,n_structure_zone,
     &          maxn_structure_zone,lmax,lmax_computed,
     &          source_type
        integer top_fluid,bot_fluid
        real*8 time_series_length,omega_imag
        integer save_velo
        complex*16 station_disp_solid(3,max_nstation)
        complex*16 station_disp_fluid(3,max_nstation)
        real*8 rmin_structure_zone(maxn_structure_zone)
        real*8 rmax_structure_zone(maxn_structure_zone)
        real*8 rho_structure_zone(4,maxn_structure_zone)
        real*8 vpv_structure_zone(4,maxn_structure_zone)
        real*8 vph_structure_zone(4,maxn_structure_zone)
        real*8 vsv_structure_zone(4,maxn_structure_zone)
        real*8 vsh_structure_zone(4,maxn_structure_zone)
        real*8 eta_structure_zone(4,maxn_structure_zone)
        complex*16 epsilon11(max_nstation),
     &            epsilon22(max_nstation),
     &            epsilon33(max_nstation),
     &            epsilon12(max_nstation),
     &            epsilon13(max_nstation),
     &            epsilon23(max_nstation)
        complex*16 pressure(max_nstation)
        complex*16 potential(max_nstation)
        complex*16 coef_c(0:lmax,-2:2,2,3)
        complex*16 coef_dcdr(0:lmax,-2:2,2,3)
        integer ir_dep_solid(max_ndep)
        real*8 r_sta_solid(max_ndep),r_sta_fluid(max_ndep)
        integer izone_idep_solid(max_ndep)
        integer izone_idep_fluid(max_ndep)
        real*8 Atop,Ctop,Ftop,Ltop,Ntop
        real*8 Abot,Cbot,Fbot,Lbot,Nbot

c other variables
        real*8 pi,km
        parameter ( pi=3.1415926535897932d0 )
        parameter ( km=1000.0 )
        character*80 outfile_disp,outfile_velo,outfile_stress,
     &               outfile_epsilon,outfile_pressure,
     &               outfile_potential,outfile_coef_cAnddcdr
        integer istation,istation_solid,istation_fluid,icomp,
     &          idep,itheta,ir,ilay_solid,ios
        integer m,max_m,min_m
        real*8 rho_sta,vpv_sta,vph_sta,vsv_sta,vsh_sta,eta_sta
        integer izone,izone1,izone2,izone_used
        real*8 A,C,F,L,N
        real*8 lamda_fluid
        complex*16 iMulOmega
        complex*16 sigma(6)
        real*8 cal_PREM_structure
        integer ncomp_foridep(2)
        integer ncomp_fluid,ncomp_solid


        !km to meter
        station_disp_solid(:,:)=station_disp_solid(:,:)*km
        station_disp_fluid(:,:)=station_disp_fluid(:,:)*km

!             iMulOmega=dcmplx(omega_imag,2.d0*pi*
!     &              dble(i_frequency)/dble(time_series_length))
             iMulOmega=dcmplx(0,2.d0*pi*
     &              dble(i_frequency)/dble(time_series_length))


!save displacement
        if(save_velo.eq.1) then
          write(outfile_disp,
     &     "('./OUTPUT_GREENS/disp_solid/freq_',i5.5)")
     &                  i_frequency
          outfile_disp=outfile_disp//"_disp"

!        open(unit=16,file=trim(outfile_disp),status='unknown',
!     &       form='formatted',iostat=ios)
          open(unit=16,file=trim(outfile_disp),status='unknown',
     &         form='unformatted',iostat=ios)
          if(ios / = 0) stop 'error saving GREENS'
          do 900 icomp=1,3
            write(16) station_disp_solid(icomp,1:n_station_solid)
 900     continue

          close(16)
!save velocity
        else if(save_velo.eq.2) then

          write(outfile_velo,
     &     "('./OUTPUT_GREENS/velo_solid/freq_',i5.5)")
     &                  i_frequency
          outfile_velo=outfile_velo//"_velo"
          open(unit=16,file=trim(outfile_velo),status='unknown',
     &         form='unformatted',iostat=ios)
          if(ios / = 0) stop 'error saving GREENS'
          do 901 icomp=1,3
             do 1001 istation=1,n_station_solid
               station_disp_solid(icomp,istation)=
     &           station_disp_solid(icomp,istation)*iMulOmega
 1001        continue
             write(16) station_disp_solid(icomp,1:n_station_solid)
  901     continue
          close(16)
        else
          stop 'save_velo should be 1 or 2'
        end if

        write(outfile_disp,
     &     "('./OUTPUT_GREENS/disp_fluid/freq_',i5.5)")
     &                  i_frequency
        outfile_disp=outfile_disp//"_disp"
        open(unit=16,file=trim(outfile_disp),status='unknown',
     &       form='unformatted',iostat=ios)
        if(ios / = 0) stop 'error saving GREENS'
        do 1110 icomp=1,3
           write(16) station_disp_fluid(icomp,1:n_station_fluid)
 1110   continue
        close(16)

        write(outfile_stress,
     &      "('./OUTPUT_GREENS/sigma/freq_',i5.5)")
     &                    i_frequency
          outfile_stress=outfile_stress//"_sigma"
!        open(unit=15,file=trim(outfile_stress),status='unknown',
!     &         form='formatted',iostat=ios)
        open(unit=15,file=trim(outfile_stress),status='unknown',
     &         form='unformatted',iostat=ios)
        if(ios / = 0) stop 'error saving GREENS'
        do 171 idep=1,ndep_solid
          izone1 = 0
          izone2 = 0
          do 150 izone=1,n_structure_zone
            if ( rmin_structure_zone(izone).le.
     &           r_sta_solid(idep)+1.e-10 )
     &        izone1 = izone
  150     continue
          do 160 izone=n_structure_zone,1,-1
            if ( rmax_structure_zone(izone).ge.
     &                r_sta_solid(idep)-1.e-10 )
     &         izone2 = izone
  160     continue

!discontinuity
          if(izone1.ne.izone2) then
             if(vsh_structure_zone(1,izone1).lt.1.e-10.and.
     &          vsh_structure_zone(1,izone2).ne.0.d0) then
                izone_used=izone2
             else
                izone_used=izone1
             end if
          else 
             izone_used=izone1
          end if
          if(vsh_structure_zone(1,izone_used)*
     &       vsv_structure_zone(1,izone_used).le.1.e-10) then
               print *,'izone',izone_used,izone1,izone2
               stop 'Error, station_solid located in fluid media!!'
          end if
          izone_used=izone_idep_solid(idep)

!          ir=ir_dep_solid(idep)
!          if(grid_mu(1,ir).gt.1.e-4) then
!             ilay_solid=1
!          else if(grid_mu(2,ir).gt.1.e-4) then
!             ilay_solid=2
!          else if(ir.eq.ngrid_r) then
! grid_mu(1,ngrid_r) and grid_mu(2,ngrid_r) is not defined or zero! 
! We use grid_mu(2,ngrid_r-1) instead.
!             ir=ir-1
!             ilay_solid=2
!          else
!             print *,'mu=',grid_mu(1,ir),grid_mu(2,ir),ir
!             stop 'Error, mu=0.0 for solid media!'
!          end if
          rho_sta 
     &      = cal_PREM_structure( r_sta_solid(idep),
     &                        rho_structure_zone(1,izone_used))
          vpv_sta
     &      = cal_PREM_structure( r_sta_solid(idep),
     &                        vpv_structure_zone(1,izone_used))
          vph_sta
     &      = cal_PREM_structure( r_sta_solid(idep),
     &                        vph_structure_zone(1,izone_used))
          vsv_sta
     &      = cal_PREM_structure( r_sta_solid(idep),
     &                        vsv_structure_zone(1,izone_used))
          vsh_sta
     &      = cal_PREM_structure( r_sta_solid(idep),
     &                        vsh_structure_zone(1,izone_used))
          eta_sta
     &      = cal_PREM_structure( r_sta_solid(idep),
     &                        eta_structure_zone(1,izone_used))


          A=rho_sta*vph_sta*vph_sta
          C=rho_sta*vpv_sta*vpv_sta
          L=rho_sta*vsv_sta*vsv_sta
          N=rho_sta*vsh_sta*vsh_sta
          F=eta_sta*(A-2.d0*L)

!Gpa to pa
          A=A*1.e9
          C=C*1.e9
          F=F*1.e9
          L=L*1.e9
          N=N*1.e9

          if(top_fluid.ne.1.and.idep.eq.1) then
            Atop=A
            Ctop=C
            Ftop=F
            Ltop=L      
            Ntop=N
          end if
          if(bot_fluid.ne.1.and.idep.eq.ndep_solid) then
            Abot=A
            Cbot=C
            Fbot=F
            Lbot=L      
            Nbot=N
          end if

          do 172 itheta=1,ntheta
            istation_solid=(idep-1)*ntheta+itheta
            sigma(1)=epsilon11(istation_solid)*C+
     &               (epsilon22(istation_solid)+
     &                epsilon33(istation_solid))*F
            sigma(2)=epsilon11(istation_solid)*F+
     &               epsilon22(istation_solid)*A+
     &               epsilon33(istation_solid)*(A-2*N)
            sigma(3)=epsilon11(istation_solid)*F+
     &               epsilon22(istation_solid)*(A-2*N)+
     &               epsilon33(istation_solid)*A
            sigma(4)=2.0*epsilon12(istation_solid)*L
            sigma(5)=2.0*epsilon13(istation_solid)*L
            sigma(6)=2.0*epsilon23(istation_solid)*N

! copy stress to array and then save in file.
            epsilon11(istation_solid)=sigma(1)
            epsilon22(istation_solid)=sigma(2)
            epsilon33(istation_solid)=sigma(3)
            epsilon12(istation_solid)=sigma(4)
            epsilon13(istation_solid)=sigma(5)
            epsilon23(istation_solid)=sigma(6)
 172      end do
 171    end do
        write(15)epsilon11(1:n_station_solid)
        write(15)epsilon22(1:n_station_solid)
        write(15)epsilon33(1:n_station_solid)
        write(15)epsilon12(1:n_station_solid)
        write(15)epsilon13(1:n_station_solid)
        write(15)epsilon23(1:n_station_solid)

        close(15)

!pressure = -kappa*div(displacent)

!Gpa to pa
        pressure(:)=pressure(:)*1e9
        write(outfile_pressure,
     &    "('./OUTPUT_GREENS/pressure/freq_',i5.5)")
     &                  i_frequency
        open(unit=15,file=trim(outfile_pressure),status='unknown',
     &       form='unformatted',iostat=ios)
        if(ios / = 0) stop 'error saving GREENS'
                write(15)pressure(1:n_station_fluid)
        close(15)


        if(i_frequency.ne.0) then
           do 181 istation_fluid=1,n_station_fluid
              potential(istation_fluid)=
     &                          pressure(istation_fluid)/iMulOmega
181        end do
c       avoid to devide zero
        else
              potential(:)=dcmplx(0.0,0.0)
        end if
        write(outfile_potential,
     &    "('./OUTPUT_GREENS/potential/freq_',i5.5)")
     &                  i_frequency
        open(unit=16,file=trim(outfile_potential),status='unknown',
     &       form='unformatted',iostat=ios)
!        open(unit=16,file=trim(outfile_potential),status='unknown',
!     &      iostat=ios)

        if(ios / = 0) stop 'error saving GREENS'
                write(16)potential(1:n_station_fluid)
!                write(16,*)potential(1:n_station_fluid)

        close(16)





!save the coefficients c and dcdr

        if(source_type.eq.1) then
          max_m=2
          min_m=-2
        else if(source_type.eq.2) then
          max_m=1
          min_m=-1
        end if
        if(top_fluid.eq.1) then
            ncomp_foridep(1)=1
        else
            ncomp_foridep(1)=3
        end if
        if(bot_fluid.eq.1) then
            ncomp_foridep(2)=1
        else
            ncomp_foridep(2)=3
        end if

        write(outfile_coef_cAnddcdr,
     &    "('./OUTPUT_GREENS/coef_cAnddcdr/freq_',i5.5)")
     &                  i_frequency
        open(unit=17,file=trim(outfile_coef_cAnddcdr),status='unknown',
     &       form='unformatted',iostat=ios)
        if(ios / = 0) stop 'error saving GREENS'
        write(17) lmax_computed
        do idep=1,2
          do icomp=1,ncomp_foridep(idep)
             do m=min_m,max_m
              write(17)coef_c(0:lmax_computed,m,idep,icomp)
              write(17)coef_dcdr(0:lmax_computed,m,idep,icomp)
             end do
          end do
        end do

        close(17)


        end subroutine


        subroutine save_par(n_frequency,rank,
     &         ndist_solid,dist_solid_table,ndist_fluid,
     &         dist_fluid_table,max_ntheta,
     &         ndep_solid,ndep_fluid,source_type,
     &         time_series_length,omega_imag,lmax_allfreq,
     &         top_fluid,bot_fluid,ir_dep_fluid,ir_dep_solid,
     &         max_ndep,maxngrid_r,grid_r,grid_rho,
     &         Atop,Ctop,Ftop,Ltop,Ntop,
     &         Abot,Cbot,Fbot,Lbot,Nbot)
        implicit none
c input/output
        integer n_frequency,rank,ndep_solid,ndep_fluid
        integer ndist_solid,ndist_fluid
        real*8 dist_fluid_table(max_ntheta),dist_solid_table(max_ntheta)
        real*8 time_series_length,omega_imag
        integer lmax_allfreq
        integer source_type
        integer top_fluid,bot_fluid
        integer max_ndep,max_ntheta,maxngrid_r
c        real*8 r_sta_solid(max_ndep),r_sta_fluid(max_ndep)
        integer ir_dep_solid(max_ndep),ir_dep_fluid(max_ndep)
        real*8 grid_rho(2,maxngrid_r),grid_r(maxngrid_r)
        real*8  Atop,Ctop,Ftop,Ltop,Ntop
        real*8  Abot,Cbot,Fbot,Lbot,Nbot

c other variables
        character*80 outfile_par
        integer idep,ios
        integer ndist,ndist_alldep_tmp
        real*8 dist_min,ddist
        real*8 pi
        parameter ( pi=3.1415926535897932d0 )
        
c        integer nl_eachpack



          if(ndep_solid.gt.0) then
             dist_min=dist_solid_table(1)
             ddist=dist_solid_table(2)-dist_solid_table(1)
             ndist=ndist_solid
          else
             dist_min=dist_fluid_table(1)
             ddist=dist_fluid_table(2)-dist_fluid_table(1)
             ndist=ndist_fluid
          end if

          dist_min=dist_min*180.0/pi
          ddist=ddist*180.0/pi

          outfile_par="./OUTPUT_GREENS/Green_Par_forConvolution"
          open(unit=18,file=trim(outfile_par),status='unknown',
     &       form='formatted',iostat=ios)
          if(ios / = 0) stop 'Error in openning file 
     &                        Green_Par_forConvolution'
          write(18,*)n_frequency
          write(18,*)0.0,0.1,ndep_fluid
          write(18,*)0.0,0.1,ndep_solid
          write(18,*)dist_min,ddist,ndist
          write(18,*)time_series_length
          write(18,*)omega_imag
          close(18)


*****************************Double couple mechanism ******************
*************************deep earth or receiver side coupling**********
         if(source_type.eq.1) then
            outfile_par="./OUTPUT_GREENS/topbot_forResave"
            open(unit=19,file=trim(outfile_par),status='unknown',
     &       form='formatted',iostat=ios)
              if(top_fluid.eq.1) then
                 write(19,*) '.true.'
              else
                 write(19,*) '.false.'
              end if
              if(bot_fluid.eq.1) then
                 write(19,*) '.true.'
              else
                 write(19,*) '.false.'
              end if
              write(19,*) 1
              write(19,*) n_frequency
c Including l=0, so it is lmax_allfreq+1
              write(19,*) lmax_allfreq+1
              write(19,*) 'nl_eachpack'
              close(19)

            outfile_par="./OUTPUT_GREENS/Partop_forSpecCoefToTime"
            open(unit=20,file=trim(outfile_par),status='unknown',
     &       form='formatted',iostat=ios)
              if(top_fluid.eq.1) then
                 write(20,*) '.true.'
              else
                 write(20,*) '.false.'
              end if
              write(20,*) source_type
              write(20,*) "ndist_topbound 1"  
              write(20,*) lmax_allfreq,'nl_eachpack',2
              write(20,*) n_frequency
              write(20,*) time_series_length
              write(20,*) omega_imag
              write(20,*) 'tbegin_cut'
              write(20,*) 'tend_cut'
              write(20,*) 'npt_eachpack'
c             write(20,*) 'f1 f2'
              close(20)

            outfile_par="./OUTPUT_GREENS/matetop_forSpecCoefToTime"
            open(unit=21,file=trim(outfile_par),status='unknown',
     &       form='formatted',iostat=ios)
            if(top_fluid.eq.1) then
              write(21,*) grid_r(ir_dep_fluid(1))
              write(21,*) grid_rho(1,ir_dep_fluid(1))
            else
              write(21,*) grid_r(ir_dep_solid(1))
              write(21,*) Atop
              write(21,*) Ctop
              write(21,*) Ftop
              write(21,*) Ltop
              write(21,*) Ntop
            end if
            close(21)

            outfile_par="./OUTPUT_GREENS/Parbot_forSpecCoefToTime"
            open(unit=22,file=trim(outfile_par),status='unknown',
     &       form='formatted',iostat=ios)
              if(bot_fluid.eq.1) then
                 write(22,*) '.true.'
              else
                 write(22,*) '.false.'
              end if
              write(22,*) source_type
              write(22,*) "ndist_botbound 1"
              write(22,*) lmax_allfreq,'nl_eachpack',2
              write(22,*) n_frequency
              write(22,*) 'nl_eachpack'
              write(22,*) time_series_length
              write(22,*) omega_imag
              write(22,*) 'tbegin_cut'
              write(22,*) 'tend_cut'
              write(22,*) 'npt_eachpack'
c             write(22,*) 'f1 f2'
              close(22)


            outfile_par="./OUTPUT_GREENS/matebot_forSpecCoefToTime"
            open(unit=23,file=trim(outfile_par),status='unknown',
     &       form='formatted',iostat=ios)
            if(bot_fluid.eq.1) then
              write(23,*) grid_r(ir_dep_fluid(ndep_fluid))
              write(23,*) grid_rho(1,ir_dep_fluid(ndep_fluid))
            else
              write(23,*) grid_r(ir_dep_solid(ndep_solid))
              write(23,*) Abot
              write(23,*) Cbot
              write(23,*) Fbot
              write(23,*) Lbot
              write(23,*) Nbot
            end if
            close(23)

           if(ndep_fluid.gt.0) then
            ndist_alldep_tmp=ndep_fluid*ndist_fluid
            outfile_par="./OUTPUT_GREENS/ParChidot_forSpectotimeEdge"
            open(unit=24,file=trim(outfile_par),status='unknown',
     &       form='formatted',iostat=ios)
            write(24,*) ndist_alldep_tmp,ndep_fluid,ndist_fluid
            write(24,*) n_frequency
            write(24,*) time_series_length
            write(24,*) "1"
            write(24,*) omega_imag
            write(24,*) 'tbegin_save'
            write(24,*) 'tend_save'
            write(24,*) 'npt_eachpack'
            close(24)


            outfile_par="./OUTPUT_GREENS/ParDispfluid_forSpectotimeEdge"
            open(unit=25,file=trim(outfile_par),status='unknown',
     &       form='formatted',iostat=ios)
            write(25,*) ndist_alldep_tmp,ndep_fluid,ndist_fluid
            write(25,*) n_frequency
            write(25,*) time_series_length
            write(25,*) '3'
            write(25,*) omega_imag
            write(25,*) 'tbegin_save'
            write(25,*) 'tend_save'
            write(25,*) 'npt_eachpack'
            close(25)
           end if

           if(ndep_solid.gt.0) then
            ndist_alldep_tmp=ndep_solid*ndist_solid
            outfile_par="./OUTPUT_GREENS/ParStress_forSpectotimeEdge"
            open(unit=26,file=trim(outfile_par),status='unknown',
     &       form='formatted',iostat=ios)
            write(26,*) ndist_alldep_tmp,ndep_solid,ndist_solid
            write(26,*) n_frequency
            write(26,*) time_series_length
            write(26,*) '6'
            write(26,*) omega_imag
            write(26,*) 'tbegin_save'
            write(26,*) 'tend_save'
            write(26,*) 'npt_eachpack'
            close(26)


            outfile_par="./OUTPUT_GREENS/ParVelosolid_forSpectotimeEdge"
            open(unit=27,file=trim(outfile_par),status='unknown',
     &       form='formatted',iostat=ios)
            write(27,*) ndist_alldep_tmp,ndep_solid,ndist_solid
            write(27,*) n_frequency
            write(27,*) time_series_length
            write(27,*) '3'
            write(27,*) omega_imag
            write(27,*) 'tbegin_save'
            write(27,*) 'tend_save'
            write(27,*) 'npt_eachpack'
            close(27)
           end if

         end if

        end subroutine save_par
