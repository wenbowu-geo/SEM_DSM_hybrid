subroutine read_parameter_SEMtoTele(Imain_SEMtoTele,nx_TopoTaper,ny_TopoTaper,nx_notopo,&
                         ny_notopo,Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh, &
                         Tele_irTop,Tele_irBot,NEX_XI,NEX_ETA,NER,NDOUBLINGS,USE_REGULAR_MESH,&
                         ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
                         CENTER_LATITUDE_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,&      
                         rotation_matrix,ner_doublings,myrank)

implicit none
include "constants.h"

integer::nx_TopoTaper,ny_TopoTaper,nx_notopo,ny_notopo
integer::Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,Tele_irTop,Tele_irBot
integer::NEX_XI,NEX_ETA,NER,NDOUBLINGS
double precision ::ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH

!WENBO
! rotation matrix from Euler angles
double precision, dimension(NDIM,NDIM) :: rotation_matrix

logical USE_REGULAR_MESH
integer::ner_doublings(2)
integer::myrank
integer::Imain_SEMtoTele

integer, external :: err_occurred


! open parameter file
 open(unit=IIN,file=MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES)) &
       //'SEMtoTele_Par_file',status='old',action='read')

 call read_value_double_precision(IIN,IGNORE_JUNK,ANGULAR_WIDTH_XI_IN_DEGREES,'SEMtoTele.ANGULAR_WIDTH_XI_IN_DEGREES')
 if(err_occurred() /= 0) stop 'Error reading SEMtoTele parameter ANGULAR_WIDTH_XI_IN_DEGREES'

 call read_value_double_precision(IIN,IGNORE_JUNK,ANGULAR_WIDTH_ETA_IN_DEGREES,'SEMtoTele.ANGULAR_WIDTH_ETA_IN_DEGREES')
 if(err_occurred() /= 0) stop 'Error reading SEMtoTele parameter ANGULAR_WIDTH_ETA_IN_DEGREES'

 call read_value_double_precision(IIN,IGNORE_JUNK,CENTER_LATITUDE_IN_DEGREES,'SEMtoTele.CENTER_LATITUDE_IN_DEGREES')
 if(err_occurred() /= 0) stop 'Error reading SEMtoTele parameter CENTER_LATITUDE_IN_DEGREES'

 call read_value_double_precision(IIN,IGNORE_JUNK,CENTER_LONGITUDE_IN_DEGREES,'SEMtoTele.CENTER_LONGITUDE_IN_DEGREES')
 if(err_occurred() /= 0) stop 'Error reading SEMtoTele parameter CENTER_LONGITUDE_IN_DEGREES'

 call read_value_double_precision(IIN,IGNORE_JUNK,GAMMA_ROTATION_AZIMUTH,'SEMtoTele.GAMMA_ROTATION_AZIMUTH')
 if(err_occurred() /= 0) stop 'Error reading SEMtoTele parameter GAMMA_ROTATION_AZIMUTH'

!WENBO
 call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)
!WENBO


 call read_value_integer(IIN,IGNORE_JUNK,nx_TopoTaper, 'SEMtoTele.nx_TopographyTaper')
 if(err_occurred() /= 0) return
 call read_value_integer(IIN,IGNORE_JUNK,ny_TopoTaper, 'SEMtoTele.ny_TopographyTaper')
 if(err_occurred() /= 0) return
 call read_value_integer(IIN,IGNORE_JUNK,nx_notopo, 'SEMtoTele.nx_Notopograpgy')
 if(err_occurred() /= 0) return
 call read_value_integer(IIN,IGNORE_JUNK,ny_notopo, 'SEMtoTele.ny_Notopograpgy')
 if(err_occurred() /= 0) return
 call read_value_integer(IIN,IGNORE_JUNK,Tele_ixLow, 'SEMtoTele.ix_LowBound')
 if(err_occurred() /= 0) return
 call read_value_integer(IIN,IGNORE_JUNK,Tele_ixHigh, 'SEMtoTele.ix_HighBound')
 if(err_occurred() /= 0) return
 call read_value_integer(IIN,IGNORE_JUNK,Tele_iyLow, 'SEMtoTele.iy_LowBound')
 if(err_occurred() /= 0) return
 call read_value_integer(IIN,IGNORE_JUNK,Tele_iyHigh, 'SEMtoTele.iy_HighBound')
 if(err_occurred() /= 0) return
 call read_value_integer(IIN,IGNORE_JUNK,Tele_irTop, 'SEMtoTele.ir_Top')
 if(err_occurred() /= 0) return
 call read_value_integer(IIN,IGNORE_JUNK,Tele_irBot, 'SEMtoTele.ir_Bound')
 if(err_occurred() /= 0) return


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


 if(Tele_ixLow.lt.0.or.Tele_ixHigh.gt.NEX_XI.or.Tele_iyLow.lt.0 &
    .or.Tele_iyHigh.gt.NEX_ETA.or.Tele_irTop.gt.NER.or.Tele_irTop.lt.0&
    .or.Tele_irBot.gt.NER.or.Tele_irBot.lt.0) then
      print *,Tele_ixLow,Tele_ixHigh,NEX_XI,Tele_iyLow,Tele_iyHigh,NEX_ETA,Tele_irBot,NER
       stop 'The required Boundary elements are not inside simulating area'
 end if


 if(Tele_ixLow.gt.Tele_ixHigh.or.Tele_iyLow.gt.Tele_iyHigh) then
       stop 'The Low value should be greater than High value'
 end if

 if(Tele_ixLow.ge.nx_notopo+1.or.Tele_iyLow.ge.ny_notopo+1.or.&
    (NEX_XI-Tele_ixHigh).ge.nx_notopo+1.or.(NEX_ETA-Tele_iyHigh).ge.ny_notopo+1) then
    write(Imain_SEMtoTele,*)  'WARNING:coupling boundaries are with topography!!!!!'
     stop 'Error:coupling boundaries are with topography!!!!!'
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
      print *,Tele_ixLow,Tele_ixHigh,NEX_XI,Tele_iyLow,Tele_iyHigh,NEX_ETA,Tele_irBot,NER
       stop 'The required Boundary elements are not inside simulating area'
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
end subroutine read_parameter_SEMtoTele
