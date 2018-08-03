subroutine add_sediment(x_eval,y_eval,z_eval,Radi,flag_media,fluid_thispoint,vp_final, &
                        vs_final,rho_final)
!use tomography, only:ORIG_LAT,SPACING_LAT,ORIG_LON,SPACING_LON,topography,R_EARTH_SURF,&
!                NTRENCH,lat_trench,lon_trench
use sediment
implicit none
 !sediment
  real(kind=CUSTOM_REAL) ::theta_surf,phi_surf
  integer ::ilat,ilon
  real(kind=CUSTOM_REAL) ::topo

!set up sediment close to trench
  real(kind=CUSTOM_REAL) ::dist_min,dist,azim,bazim
  real(kind=CUSTOM_REAL) ::dist_Ocbot,dist_Ocbot_square,dist_Ocbot_cube
  real(kind=CUSTOM_REAL) ::sedi_thick_thispoint
  real(kind=CUSTOM_REAL) ::depth_slabtop
  integer ::itrench,itrench_find,sediment_units,vp_grd_accreprism,&
            taper_prism,vp_prism


!Start-         South America - sediment
!The maximum sediment thickness is less than 4000 meters
!Radi might be changed, so recalculate it.
              Radi=dsqrt(x_eval**2+y_eval**2+z_eval**2)
              theta_surf=dacos(z_eval/Radi)
              phi_surf=dacos(dble(x_eval)/(dble(Radi)*dsin(theta_surf))*(1-1.e-7))
              if(y_eval<-1.e-14) then
                    phi_surf=-phi_surf
              end if
              theta_surf=theta_surf*180.0/3.1415926
              phi_surf=phi_surf*180.0/3.1415926
              ilat=(90.0-theta_surf-ORIG_LAT)/SPACING_LAT
              ilon=(phi_surf-ORIG_LON)/SPACING_LON
              topo=topography(ilat,ilon)

    if(Radi-R_EARTH_SURF-topo.gt.-15000.0) then
              dist_min=180.0

              do itrench=1,NTRENCH
                     call distazbaz(90.0-theta_surf,phi_surf,lat_trench(itrench),&
                       lon_trench(itrench),dist,azim,bazim)
                       if(dist<dist_min) then
                          dist_min=dist
                          itrench_find=itrench
                       end if
              end do



!              if(myrank.eq.55) print
!              *,'dist=',dist_min,90.0-theta_surf,phi_surf,lat_trench(itrench_find),&
!                lon_trench(itrench_find)
!sediment layer. Vs, Vp and rho relationship is refered to Thomas M. Brocher,
!BSSA(2010)
!              if(Radi-R_EARTH_SURF-topo.gt.-3000.0.and.&
!                 topo.lt.0.and.fluid_thispoint.eq.0&
!                 .and.phi_surf.ge.lon_trench(itrench_find)) then
!                   vs_final=1000.0-(Radi-R_EARTH_SURF-topo)/3000.0*2000.0
!                   if(vs_final.lt.1000.0) vs_final=1000.0
!                   vs_final=vs_final/1000.0
!                   vp_final=1.16*vs_final+1.36
!                   if(vp_final.gt.4.0) then
!                     vp_final=4.0
!                     vs_final=(vp_final-1.36)/1.16
!                   end if
!
!                   rho_final=1.6612*vp_final-0.4721*vp_final**2+0.0671*vp_final**3&
!                             -0.0043*vp_final**4+0.000106*vp_final**5
!                   vs_final=vs_final*1000.0
!                   vp_final=vp_final*1000.0
!                   rho_final=rho_final*1000.0
!              end if
!Marine  sediment
!Edwin L. Hmilton, JASA, 1979 and JSP, 1976
!              if(Radi-R_EARTH_SURF-topo.gt.-2300.0.and.&
!                 topo.lt.0.and.fluid_thispoint.eq.0&
!                 .and.dist_min.le.0.3) then
!              if(Radi-R_EARTH_SURF-topo.gt.-1500.0*(0.27-dist_min)/0.27.and.&
!                 topo.lt.0.and.fluid_thispoint.eq.0&
!                 .and.dist_min.le.0.27.and.90.0-theta_surf.lt.-33.0.and.&
!                 phi_surf.ge.lon_trench(itrench_find)) then

!Trench sediment
              if(fluid_thispoint.eq.0.and.dist_min.le.0.224) then
                  sedi_thick_thispoint=2300*dist_min/0.224
                  sediment_units=1
              else
if(fluid_thispoint.eq.0.and.dist_min.ge.0.224.and.topo.lt.-500) then
!Sediemnt on shelft
                  if(dist_min.ge.0.5396) then
                    sedi_thick_thispoint=1000.0
                    sediment_units=1
!                  sedi_thick_thispoint=1500.0
!Accretion prism
                  else
!1000meters layer 1. Layer 2 has linearly increased thickness, from 0 at
!distance 0.224dge to 3900meters at distance 0.5396deg.
                    depth_slabtop=7431.0+(dist_min-0.224)*111.2*1000.0*tan(10.0/180.0*3.1415926)
                    sedi_thick_thispoint=topo+depth_slabtop
!                    if(sedi_thick_thispoint>13000.0) print
!                    *,'check_sedi_thickness',dist_min,sedi_thick_thispoint
                    vp_grd_accreprism=(5800.0-5000.0)/10000.0
                    sediment_units=2

                    sedi_thick_thispoint=1000.0
                    sediment_units=1
                  end if
              else
                  sedi_thick_thispoint=0.0
              end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!No sediment
!              sedi_thick_thispoint=0.0


              if(Radi-R_EARTH_SURF-topo.gt.-sedi_thick_thispoint.and.&
                 sedi_thick_thispoint.gt.1.0.and.topo.lt.0.and.fluid_thispoint.eq.0&
                 .and.90.0-theta_surf.lt.-33.0.and.&
                 phi_surf.ge.lon_trench(itrench_find).and.flag_media.eq.2) then

                   dist_Ocbot=-(Radi-R_EARTH_SURF-topo)/1000.0

!                   if(sedi_thick_thispoint>13000.0.and.dist_Ocbot>10000.0)
!                   print *,'check_sedi_thickness',dist_min,&
!                        sediment_units,dist_Ocbot*1000.0,sedi_thick_thispoint

                   if(dist_Ocbot.lt.0) dist_Ocbot=0.0

!top sediment layer
                   if(sediment_units.eq.1.or.(sediment_units.eq.2.and.dist_Ocbot*1000.0.lt.1000)) then
                      dist_Ocbot_square=dist_Ocbot*dist_Ocbot
                      dist_Ocbot_cube=dist_Ocbot_square*dist_Ocbot
                      if(dist_Ocbot.lt.0) print *,'dist_Ocbot',dist_Ocbot,dist_Ocbot_square,dist_Ocbot_cube,&
                         topo,90.0-theta_surf,phi_surf
                      if(dist_Ocbot<1.0) then
                        vp_final=1.511+dist_Ocbot*1.304-0.741*dist_Ocbot_square+&
                            0.257*dist_Ocbot_cube
                        if(vp_final.gt.4.0) then
                          vp_final=4.0
                        end if
                        vs_final=0.78*vp_final-0.962
                        rho_final=1.53+1.395*dist_Ocbot-0.617*dist_Ocbot_square
                      else
                        vp_final=2.331+0.593*(dist_Ocbot-1.0)
                        if(vp_final.gt.4.0) then
                          vp_final=4.0
                        end if
                        vs_final=0.78*vp_final-0.962
                        rho_final=2.308+0.161*(dist_Ocbot-1)
                      end if
                      vs_final=vs_final*1000.0
                      vp_final=vp_final*1000.0
                      rho_final=rho_final*1000.0

                      if(vs_final.lt.500.0) vs_final=500.0
!the second layer of accretion prism
                    else
                      if(dist_min<0.4496) then
                        taper_prism=1.0
                      else
                        taper_prism=(0.5396-dist_min)/(0.5396-0.4496)
                      end if
                      vp_prism=5000.0+(dist_Ocbot*1000.0-1000.0)*vp_grd_accreprism
!                      if(Radi-R_EARTH_SURF<-13000.0) &
                         print *, 'error check_prism',Radi-R_EARTH_SURF,vp_prism
                      vp_final=vp_prism*taper_prism+5800.0*(1-taper_prism)
                      if(vp_final>5800.0) vp_final=5800.0
                      vs_final=vp_final/1.9
                      rho_final=2600.0
                    end if

              end if
      end if !Radi-R_EARTH_SURF-topo.gt.-15000.0
!!!!!!END-  South America
!sediment********************************************************************
  end if !sediment incorporated or not


end subroutine add_sediment
