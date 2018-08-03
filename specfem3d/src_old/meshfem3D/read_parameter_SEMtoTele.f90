subroutine read_parameter_SEMtoTele(myrank)

!include "constants.h"

use constants
use SEMtoTele_MeshPar
!use meshfem3D_par, only:Z_DEPTH_BLOCK,NEX_XI,NEX_ETA,NER,&
!        NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NDOUBLINGS,&
!        USE_REGULAR_MESH,ner_doublings,iproc_xi_current,iproc_eta_current
implicit none


integer myrank
!integer::nx_TopoTaper,ny_TopoTaper,nx_notopo,ny_notopo
!integer::Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,Tele_irTop,Tele_irBot
!integer::NEX_XI,NEX_ETA,NER,NDOUBLINGS
!logical USE_REGULAR_MESH
!integer::ner_doublings(2)
!integer::myrank
!integer::Imain_SEMtoTele
!double precision::x_temp(NGNOD_EIGHT_CORNERS),y_temp(NGNOD_EIGHT_CORNERS),&
!                   z_temp(NGNOD_EIGHT_CORNERS)
!double precision ::radi,theta,phi,lat_corner,lon_corner
integer :: ier


 if(myrank == 0 .and. Imain_SEMtoTele /= ISTANDARD_OUTPUT) &
       open(unit=Imain_SEMtoTele,file=trim(OUTPUT_FILES)//'/output_SEMtoTele.txt',status='unknown')
 if(myrank.eq.0) then
      write(Imain_SEMtoTele,*)
      write(Imain_SEMtoTele,*) '******************************************'
      write(Imain_SEMtoTele,*) '*** INTERFACING WITH SEM (WENBO WU)***'
      write(Imain_SEMtoTele,*) '******************************************'
      write(Imain_SEMtoTele,*)
      write(Imain_SEMtoTele,*) 'Reading the interfacing boundary information from file ',&
                MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES))//'SEMtoTele_Par_file'
 end if


! open parameter file
 open(unit=IIN,file=MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES)) &
       //'SEMtoTele_Par_file',status='old',action='read')
 call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,ANGULAR_WIDTH_XI_IN_DEGREES,'SEMtoTele.ANGULAR_WIDTH_XI_IN_DEGREES',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter ANGULAR_WIDTH_XI_IN_DEGREES'

 call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,ANGULAR_WIDTH_ETA_IN_DEGREES,'SEMtoTele.ANGULAR_WIDTH_ETA_IN_DEGREES',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter ANGULAR_WIDTH_ETA_IN_DEGREES'

 call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,CENTER_LATITUDE_IN_DEGREES,'SEMtoTele.CENTER_LATITUDE_IN_DEGREES',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter CENTER_LATITUDE_IN_DEGREES'

 call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,CENTER_LONGITUDE_IN_DEGREES,'SEMtoTele.CENTER_LONGITUDE_IN_DEGREES',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter CENTER_LONGITUDE_IN_DEGREES'

 call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,GAMMA_ROTATION_AZIMUTH,'SEMtoTele.GAMMA_ROTATION_AZIMUTH',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter GAMMA_ROTATION_AZIMUTH'

! call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,DEPTH_BLOCK_KM_tmp,'SEMtoTele.DEPTH_BLOCK_KM',ier)
! if(ier /= 0) stop 'Error reading SEMtoTele paramete DEPTH_BLOCK_KM'

! Z_DEPTH_BLOCK = - dabs(DEPTH_BLOCK_KM_tmp) * 1000.d0


 call read_value_integer_mesh(IIN,IGNORE_JUNK,NX_TOPOTAPER, 'SEMtoTele.nx_TopographyTaper',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter NX_TOPOTAPER'

 call read_value_integer_mesh(IIN,IGNORE_JUNK,NY_TOPOTAPER, 'SEMtoTele.ny_TopographyTaper',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter NY_TOPOTAPER'

 call read_value_integer_mesh(IIN,IGNORE_JUNK,NX_NOTOPO, 'SEMtoTele.nx_Notopograpgy',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter NX_NOTOPO'

 call read_value_integer_mesh(IIN,IGNORE_JUNK,NY_NOTOPO, 'SEMtoTele.ny_Notopograpgy',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter NY_NOTOPO'

 call read_value_integer_mesh(IIN,IGNORE_JUNK,TELE_IXLOW, 'SEMtoTele.ix_LowBound',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IXLOW'

 call read_value_integer_mesh(IIN,IGNORE_JUNK,TELE_IXHIGH, 'SEMtoTele.ix_HighBound',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IXHIGH'

 call read_value_integer_mesh(IIN,IGNORE_JUNK,TELE_IYLOW, 'SEMtoTele.iy_LowBound',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IYLOW'

 call read_value_integer_mesh(IIN,IGNORE_JUNK,TELE_IYHIGH, 'SEMtoTele.iy_HighBound',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IYHIGH'

 call read_value_integer_mesh(IIN,IGNORE_JUNK,TELE_IRTOP, 'SEMtoTele.ir_Top',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IRTOP'

 call read_value_integer_mesh(IIN,IGNORE_JUNK,TELE_IRBOT, 'SEMtoTele.ir_Bound',ier)
 if(ier /= 0) stop 'Error reading SEMtoTele parameter TELE_IRBOT'


! close parameter file
 close(IIN)
 if(myrank.eq.0) then
     write(Imain_SEMtoTele,*) nx_TopoTaper,' x-elements are within topography taper ranges'
     write(Imain_SEMtoTele,*) nx_TopoTaper,' y-elements are within topography taper ranges'
     write(Imain_SEMtoTele,*) 'The lower-x boundary elements is xi=',Tele_ixLow,'elements'
     write(Imain_SEMtoTele,*) 'The higher-x boundary elements is  xi=',Tele_ixHigh,'elements'
     write(Imain_SEMtoTele,*) 'The lower-y boundary elements is yi=',Tele_iyLow,'elements'
     write(Imain_SEMtoTele,*) 'The higher-y boundary elements is xi=',Tele_iyHigh,'elements'
     write(Imain_SEMtoTele,*) 'The top boundary elements is ri=',Tele_irTop,'elements'
     write(Imain_SEMtoTele,*) 'The bottom boundary elements is ri=',Tele_irBot,'elements'
     write(Imain_SEMtoTele,*)
 end if
! print *,'read_TelePar',nx_TopoTaper,nx_TopoTaper,Tele_ixLow,Tele_ixHigh
end subroutine read_parameter_SEMtoTele


!**************************************************
subroutine  prepare_variables_SEMtoTele(myrank)
use constants
use SEMtoTele_MeshPar
use meshfem3D_par, only:Z_DEPTH_BLOCK,NEX_XI,NEX_ETA,NER,&
        NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NDOUBLINGS,&
        USE_REGULAR_MESH,ner_doublings,iproc_xi_current,iproc_eta_current
implicit none

integer myrank
!integer::nx_TopoTaper,ny_TopoTaper,nx_notopo,ny_notopo
!integer::Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,Tele_irTop,Tele_irBot
!integer::NEX_XI,NEX_ETA,NER,NDOUBLINGS
!logical USE_REGULAR_MESH
!integer::ner_doublings(2)
!integer::myrank
!integer::Imain_SEMtoTele
double precision::x_temp(NGNOD_EIGHT_CORNERS),y_temp(NGNOD_EIGHT_CORNERS),&
                   z_temp(NGNOD_EIGHT_CORNERS)
double precision ::radi,theta,phi,lat_corner,lon_corner
integer :: ier,icorner



 if(Tele_ixLow.lt.0.or.Tele_ixHigh.gt.NEX_XI.or.Tele_iyLow.lt.0 &
    .or.Tele_iyHigh.gt.NEX_ETA.or.Tele_irTop.gt.NER.or.Tele_irTop.lt.0&
    .or.Tele_irBot.gt.NER.or.Tele_irBot.lt.0) then
      print *,'check_TelePar_input',Tele_ixLow,Tele_ixHigh,NEX_XI,Tele_iyLow,Tele_iyHigh,NEX_ETA,&
        Tele_irBot,Tele_irTop,NER
      call exit_MPI(myrank,'The required Boundary elements are not inside simulating area')      
 end if


 if(Tele_ixLow.gt.Tele_ixHigh.or.Tele_iyLow.gt.Tele_iyHigh) then
      call exit_MPI(myrank,'The Low value should be greater than High value')
 end if

 if(Tele_ixLow.ge.nx_notopo+1.or.Tele_iyLow.ge.ny_notopo+1.or.&
    (NEX_XI-Tele_ixHigh).ge.nx_notopo+1.or.(NEX_ETA-Tele_iyHigh).ge.ny_notopo+1) then
    print *, "coup-info",Tele_ixLow,nx_notopo,Tele_iyLow,ny_notopo,NEX_XI,NEX_ETA,Tele_ixHigh,Tele_iyHigh
    write(Imain_SEMtoTele,*)  'WARNING:coupling boundaries are with topography!!!!!'
    call exit_MPI(myrank,'Error:coupling boundaries are with topography!!!!!')
 end if
 if(dabs(R_TOP_BOUND-R_EARTH_SURF).lt.TINYVAL.and.Tele_irTop.ne.NER) then
    call exit_MPI(myrank,'Error,the top boundary is free surface, Tele_irTop should &
                  be equal to NER!!!')
 end if 


 if(NDOUBLINGS.eq.1.and.(.not.USE_REGULAR_MESH)) then
!       if(modulo(Tele_ixLow-1,4).ne.0.or.modulo(Tele_ixHigh,4).ne.0.or. &
!          modulo(Tele_iyLow-1,4).ne.0.or.modulo(Tele_iyHigh,4).ne.0) then
!            call exit_MPI(myrank,'for 1 doubling layer,Tele_ixLow,Tele_ixLow,Tele_iyLow, &
!                        Tele_iyHigh should be divided evenly by 4')
!       end if
       if(modulo(Tele_ixLow-1,4).ne.0) Tele_ixLow=((Tele_ixLow-1)/4)*4+1
       if(modulo(Tele_ixHigh,4).ne.0)  Tele_ixHigh=((Tele_ixHigh-1)/4)*4
       if(modulo(Tele_iyLow-1,4).ne.0) Tele_iyLow=((Tele_iyLow-1)/4)*4+1
       if(modulo(Tele_iyHigh,4).ne.0)  Tele_iyHigh=((Tele_iyHigh-1)/4)*4
       if(Tele_irTop.eq.ner_doublings(1)) Tele_irTop=ner_doublings(1)-1
       if(Tele_irBot.eq.ner_doublings(1)) Tele_irBot=ner_doublings(1)-1
 else if(NDOUBLINGS.eq.2.and.(.not.USE_REGULAR_MESH)) then
!       if(modulo(Tele_ixLow-1,8).ne.0.or.modulo(Tele_ixHigh,8).ne.0.or. &
!          modulo(Tele_iyLow-1,8).ne.0.or.modulo(Tele_iyHigh,8).ne.0) then
!          call exit_MPI(myrank,'for 1 doubling layer,Tele_ixLow,Tele_ixLow,Tele_iyLow, &
!                        Tele_iyHigh should be divided evenly by 8')
!       end if

       if(modulo(Tele_ixLow-1,8).ne.0) Tele_ixLow=((Tele_ixLow-1)/8)*8+1
       if(modulo(Tele_ixHigh,8).ne.0)  Tele_ixHigh=((Tele_ixHigh-1)/8)*8
       if(modulo(Tele_iyLow-1,8).ne.0) Tele_iyLow=((Tele_iyLow-1)/8)*8+1
       if(modulo(Tele_iyHigh,8).ne.0)  Tele_iyHigh=((Tele_iyHigh-1)/8)*8
       if(Tele_irTop.eq.ner_doublings(1)) then
          Tele_irTop=ner_doublings(1)-1
       else if(Tele_irTop.eq.ner_doublings(2)) then
          Tele_irTop=ner_doublings(2)-1
       end if
       if(Tele_irBot.eq.ner_doublings(1)) then
          Tele_irBot=ner_doublings(1)-1
       else if(Tele_irBot.eq.ner_doublings(2)) then
          Tele_irBot=ner_doublings(2)-1
       end if

 end if

   
 if(Tele_ixLow.lt.0.or.Tele_ixHigh.gt.NEX_XI.or.Tele_iyLow.lt.0 &
    .or.Tele_iyHigh.gt.NEX_ETA.or.Tele_irTop.gt.NER.or.Tele_irTop.lt.0&
    .or.Tele_irBot.gt.NER.or.Tele_irBot.lt.0) then
      print *,'check_TelePar_input1',Tele_ixLow,Tele_ixHigh,NEX_XI,Tele_iyLow,Tele_iyHigh,&
        NEX_ETA,Tele_irTop,Tele_irBot,NER
      call exit_MPI(myrank,'The required Boundary elements are not inside simulating area')
 end if


 if(myrank.eq.0) then
    write (Imain_SEMtoTele,*) 'There ',NDOUBLINGS,'doubling layers'
    if(NDOUBLINGS.eq.1)  then
       write (Imain_SEMtoTele,*) 'The number of the xi-low and yi-low boundary element should be 4*multiple+1'
       write (Imain_SEMtoTele,*) 'The number of the xi-high and yi-high boundary element should be 4*multiple'
       write (Imain_SEMtoTele,*) 'if ri-dow=the upper element layer of the doubling layer(including 2  &
                                 element layers),make it to equal to the lower'
    else if(NDOUBLINGS.eq.2) then
       write (Imain_SEMtoTele,*) 'The number of the xi-low and yi-low boundary element should be 8*multiple+1'
       write (Imain_SEMtoTele,*) 'The number of the xi-high and yi-high boundary element should be 8*multiple'
       write (Imain_SEMtoTele,*) 'if ri-dow=the upper element layer of the doubling layer(including 2  &
                                  element layers),make it to equal to the lower'

    end if
 
    write(Imain_SEMtoTele,*)  'The new bounddary elements is xi-low=',Tele_ixLow
    write(Imain_SEMtoTele,*)  '                              xi-high=',Tele_ixHigh
    write(Imain_SEMtoTele,*)  '                              iy-Low=',Tele_iyLow
    write(Imain_SEMtoTele,*)  '                              iy-high=',Tele_iyHigh
    write(Imain_SEMtoTele,*)  '                              ri-up=',Tele_irTop

    write(Imain_SEMtoTele,*)  '                              ri-down=',Tele_irBot
    write(Imain_SEMtoTele,*)
    write(Imain_SEMtoTele,*)
    write(Imain_SEMtoTele,*)
    write(Imain_SEMtoTele,*)
 end if

 if(myrank.eq.0)  write(Imain_SEMtoTele,*) 'counting the number of the boundary elements'

 call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

 call compute_parameter_SEMtoTele(Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,&
               Tele_irTop,Tele_irBot,TeleEle_nxLow,&
               TeleEle_nxHigh,TeleEle_nyLow,TeleEle_nyHigh,TeleEle_nrTop,&
               TeleEle_nrBot,iproc_xi_current,iproc_eta_current,&
               NEX_XI,NEX_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER,&
               USE_REGULAR_MESH,NDOUBLINGS,ner_doublings,DEBUG_COUPLING,myrank)

 
 if(myrank.eq.0) then
  x_temp(1)=0.0
  y_temp(1)=0.0
  z_temp(1)=0.0

  x_temp(2)=ANGULAR_WIDTH_XI_IN_DEGREES
  y_temp(2)=0.0
  z_temp(2)=0.0

  x_temp(3)=ANGULAR_WIDTH_XI_IN_DEGREES
  y_temp(3)=ANGULAR_WIDTH_ETA_IN_DEGREES
  z_temp(3)=0.0

  x_temp(4)=0.0
  y_temp(4)=ANGULAR_WIDTH_ETA_IN_DEGREES
  z_temp(4)=0.0

  x_temp(5)=0.0
  y_temp(5)=0.0
  z_temp(5)=0.0

  x_temp(6)=ANGULAR_WIDTH_XI_IN_DEGREES
  y_temp(6)=0.0
  z_temp(6)=0.0

  x_temp(7)=ANGULAR_WIDTH_XI_IN_DEGREES
  y_temp(7)=ANGULAR_WIDTH_ETA_IN_DEGREES
  z_temp(7)=0.0

  x_temp(8)=0.0
  y_temp(8)=ANGULAR_WIDTH_ETA_IN_DEGREES
  z_temp(8)=0.0


  !compute the four corners of simulated region. 
  call cubed_geo(x_temp,y_temp,z_temp,ANGULAR_WIDTH_XI_IN_DEGREES, &
           ANGULAR_WIDTH_ETA_IN_DEGREES,rotation_matrix)

  write(Imain_SEMtoTele,*) 'this simulated cubed sphere region is:'
  do icorner=1,4
    radi=dsqrt(x_temp(icorner)**2+y_temp(icorner)**2+z_temp(icorner)**2)
    theta=dacos(z_temp(icorner)/radi)
    if(y_temp(icorner)<0.0) then
       phi=2*PI-dacos(x_temp(icorner)*(1.d0-1.e-10)/(radi*dsin(theta)))
    else
       phi=dacos(x_temp(icorner)*(1.d0-1.e-10)/(radi*dsin(theta)))
    end if
    lat_corner= 180/PI*(PI/2-theta)
    lon_corner=phi/PI*180.0
    write(Imain_SEMtoTele,*) 'corner',icorner,'  lat=',lat_corner,' lon=',lon_corner
   end do

  write(Imain_SEMtoTele,*)
  write(Imain_SEMtoTele,*)
  write(Imain_SEMtoTele,*) 'the number of the boundary elements in proc0:'
  write(Imain_SEMtoTele,*) 'number of the x-low=',TeleEle_nxLow
  write(Imain_SEMtoTele,*) 'number of the x-high=',TeleEle_nxHigh
  write(Imain_SEMtoTele,*) 'number of the y-low=',TeleEle_nyLow
  write(Imain_SEMtoTele,*) 'number of the y-high=',TeleEle_nyHigh
  write(Imain_SEMtoTele,*) 'number of the r-top=',TeleEle_nrTop
  write(Imain_SEMtoTele,*) 'number of the r-down=',TeleEle_nrBot
 end if

 allocate(TeleEle_xLow(TeleEle_nxLow))
 allocate(TeleEle_xLowReg(TeleEle_nxLow))
 allocate(TeleEle_xHigh(TeleEle_nxHigh))
 allocate(TeleEle_xHighReg(TeleEle_nxHigh))
 allocate(TeleEle_yLow(TeleEle_nyLow))
 allocate(TeleEle_yLowReg(TeleEle_nyLow))
 allocate(TeleEle_yHigh(TeleEle_nyHigh))
 allocate(TeleEle_yHighReg(TeleEle_nyHigh))
 allocate(TeleEle_rTop(TeleEle_nrTop))
 allocate(TeleEle_rTopReg(TeleEle_nrTop))
 allocate(TeleEle_rBot(TeleEle_nrBot))
 allocate(TeleEle_rBotReg(TeleEle_nrBot))
end subroutine prepare_variables_SEMtoTele
