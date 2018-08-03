!gathering the boundary element the the teleseis interfacing.
!WENBO
subroutine getEle_TeleBound(Tele_ixLow,Tele_ixHigh,Tele_iyLow, &
             Tele_iyHigh,Tele_irTop,Tele_irBot,TeleEle_nxLow,TeleEle_ixLow,TeleEle_xLow,TeleEle_xLowReg, &
             TeleEle_nxHigh,TeleEle_ixHigh,TeleEle_xHigh,TeleEle_xHighReg,TeleEle_nyLow, &
             TeleEle_iyLow,TeleEle_yLow,TeleEle_yLowReg,TeleEle_nyHigh,TeleEle_iyHigh,TeleEle_yHigh, &
             TeleEle_yHighReg,TeleEle_nrTop,TeleEle_irTop,TeleEle_rTop,TeleEle_rTopReg,&
             TeleEle_nrBot,TeleEle_irBot,TeleEle_rBot,TeleEle_rBotReg,ixele, &
             iyele,irele,isubregion,ispec_superbrick,ispec,iproc_xi,iproc_eta, &
             NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER,USE_REGULAR_MESH,NDOUBLINGS,myrank)

 implicit none

 logical USE_REGULAR_MESH
 integer NDOUBLINGS
 integer myrank
 integer ::Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,Tele_irTop,Tele_irBot
 integer ::Tele_ixLowSub,Tele_ixHighSub,Tele_iyLowSub,Tele_iyHighSub,Tele_irSub_Top,Tele_irSub_Bot

 integer ::TeleEle_nxLow,TeleEle_nxHigh,TeleEle_nyLow,TeleEle_nyHigh,TeleEle_nrTop,TeleEle_nrBot
 integer,dimension(TeleEle_nxLow)  ::TeleEle_xLow,TeleEle_xLowReg
 integer,dimension(TeleEle_nxHigh) ::TeleEle_xHigh,TeleEle_xHighReg
 integer,dimension(TeleEle_nyLow)  ::TeleEle_yLow,TeleEle_yLowReg
 integer,dimension(TeleEle_nyHigh) ::TeleEle_yHigh,TeleEle_yHighReg
 integer,dimension(TeleEle_nrTop)  ::TeleEle_rTop,TeleEle_rTopReg
 integer,dimension(TeleEle_nrBot)  ::TeleEle_rBot,TeleEle_rBotReg

  
 integer ::TeleEle_ixLow,TeleEle_ixHigh,TeleEle_iyLow,TeleEle_iyHigh,TeleEle_irTop,TeleEle_irBot
 logical ::ix_inBound,ix_onLowBound,ix_onHighBound,iy_inBound,iy_onLowBound,iy_onHighBound
 logical ::ir_inBound,ir_onTopBound,ir_onBotBound

 integer ::ixele,iyele,irele,irele_true
 integer ::isubregion,ispec_superbrick,ispec,iproc_xi,iproc_eta,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER

 integer ::int_tmp

!The below two lines are used to avoid possible error reports during compling
!the code. NER and myrank are not used now, that occationally causes error 
!reports. But they might be usful in the future, so we keep them here.
 int_tmp=NER
 int_tmp=myrank

 if(USE_REGULAR_MESH) then
   Tele_ixLowSub=(Tele_ixLow-(iproc_xi)*NEX_PER_PROC_XI)
   Tele_ixHighSub=(Tele_ixHigh-(iproc_xi)*NEX_PER_PROC_XI)
   Tele_iyLowSub=(Tele_iyLow-(iproc_eta)*NEX_PER_PROC_ETA)
   Tele_iyHighSub=(Tele_iyHigh-(iproc_eta)*NEX_PER_PROC_ETA)
   Tele_irSub_Top=Tele_irTop
   Tele_irSub_Bot=Tele_irBot
 else 
   if(NDOUBLINGS == 1) then
     select case (isubregion)
     case(1)
       Tele_ixLowSub=(Tele_ixLow-1-(iproc_xi)*NEX_PER_PROC_XI)/2+1
       Tele_ixHighSub=(Tele_ixHigh-(iproc_xi)*NEX_PER_PROC_XI)/2
       Tele_iyLowSub=(Tele_iyLow-1-(iproc_eta)*NEX_PER_PROC_ETA)/2+1
       Tele_iyHighSub=(Tele_iyHigh-(iproc_eta)*NEX_PER_PROC_ETA)/2
       Tele_irSub_Top=Tele_irTop
       Tele_irSub_Bot=Tele_irBot
     case(2)
       Tele_ixLowSub=(Tele_ixLow-1-(iproc_xi)*NEX_PER_PROC_XI)/4+1
       Tele_ixHighSub=(Tele_ixHigh-(iproc_xi)*NEX_PER_PROC_XI)/4
       Tele_iyLowSub=(Tele_iyLow-1-(iproc_eta)*NEX_PER_PROC_ETA)/4+1
       Tele_iyHighSub=(Tele_iyHigh-(iproc_eta)*NEX_PER_PROC_ETA)/4
       Tele_irSub_Top=Tele_irTop
       Tele_irSub_Bot=Tele_irBot
     case(3)
       Tele_ixLowSub=(Tele_ixLow-(iproc_xi)*NEX_PER_PROC_XI)
       Tele_ixHighSub=(Tele_ixHigh-(iproc_xi)*NEX_PER_PROC_XI)
       Tele_iyLowSub=(Tele_iyLow-(iproc_eta)*NEX_PER_PROC_ETA)
       Tele_iyHighSub=(Tele_iyHigh-(iproc_eta)*NEX_PER_PROC_ETA)
       Tele_irSub_Top=Tele_irTop
       Tele_irSub_Bot=Tele_irBot
     case default
       stop 'Wrong number of subregions'
     end select
   else if(NDOUBLINGS == 2) then
     select case (isubregion)
     case(1)
        Tele_ixLowSub=(Tele_ixLow-1-(iproc_xi)*NEX_PER_PROC_XI)/4+1
        Tele_ixHighSub=(Tele_ixHigh-(iproc_xi)*NEX_PER_PROC_XI)/4
        Tele_iyLowSub=(Tele_iyLow-1-(iproc_eta)*NEX_PER_PROC_ETA)/4+1
        Tele_iyHighSub=(Tele_iyHigh-(iproc_eta)*NEX_PER_PROC_ETA)/4
        Tele_irSub_Top=Tele_irTop
        Tele_irSub_Bot=Tele_irBot
     case(2)
        Tele_ixLowSub=(Tele_ixLow-1-(iproc_xi)*NEX_PER_PROC_XI)/8+1
        Tele_ixHighSub=(Tele_ixHigh-(iproc_xi)*NEX_PER_PROC_XI)/8
        Tele_iyLowSub=(Tele_iyLow-1-(iproc_eta)*NEX_PER_PROC_ETA)/8+1
        Tele_iyHighSub=(Tele_iyHigh-(iproc_eta)*NEX_PER_PROC_ETA)/8
        Tele_irSub_Top=Tele_irTop
        Tele_irSub_Bot=Tele_irBot
     case(3)
        Tele_ixLowSub=(Tele_ixLow-1-(iproc_xi)*NEX_PER_PROC_XI)/2+1
        Tele_ixHighSub=(Tele_ixHigh-(iproc_xi)*NEX_PER_PROC_XI)/2
        Tele_iyLowSub=(Tele_iyLow-1-(iproc_eta)*NEX_PER_PROC_ETA)/2+1
        Tele_iyHighSub=(Tele_iyHigh-(iproc_eta)*NEX_PER_PROC_ETA)/2
        Tele_irSub_Top=Tele_irTop
        Tele_irSub_Bot=Tele_irBot
      case(4)
        Tele_ixLowSub=(Tele_ixLow-1-(iproc_xi)*NEX_PER_PROC_XI)/4+1
        Tele_ixHighSub=(Tele_ixHigh-(iproc_xi)*NEX_PER_PROC_XI)/4
        Tele_iyLowSub=(Tele_iyLow-1-(iproc_eta)*NEX_PER_PROC_ETA)/4+1
        Tele_iyHighSub=(Tele_iyHigh-(iproc_eta)*NEX_PER_PROC_ETA)/4
        Tele_irSub_Top=Tele_irTop
        Tele_irSub_Bot=Tele_irBot
      case(5)
        Tele_ixLowSub=Tele_ixLow-(iproc_xi)*NEX_PER_PROC_XI
        Tele_ixHighSub=Tele_ixHigh-(iproc_xi)*NEX_PER_PROC_XI
        Tele_iyLowSub=Tele_iyLow-(iproc_eta)*NEX_PER_PROC_ETA
        Tele_iyHighSub=Tele_iyHigh-(iproc_eta)*NEX_PER_PROC_ETA
        Tele_irSub_Top=Tele_irTop
        Tele_irSub_Bot=Tele_irBot
      case default
        stop 'Wrong number of subregions'
      end select
   else
      stop 'Wrong number of doublings'
   end if
 end if
  
!debug
! print *,'start getEle',TeleEle_ixLow,TeleEle_ixHigh,TeleEle_iyLow,TeleEle_iyHigh,TeleEle_irTop, &
!         'TeleEle_nxLow',TeleEle_nxLow,TeleEle_nxHigh,TeleEle_nyLow,TeleEle_nyHigh,TeleEle_nrTop,&
!         'Tele_ixLowSub',(Tele_ixLow-(iproc_xi)*NEX_PER_PROC_XI)/2+1,Tele_ixLow,(iproc_xi),&
!        NEX_PER_PROC_XI,Tele_ixLowSub,Tele_ixHighSub,Tele_iyLowSub,Tele_iyHighSub,Tele_irSub_Top,&
!        Tele_irSub_Bot,isubregion,&
!        'Tele_ixLow',Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,Tele_irTop,Tele_irBot,myrank
 ix_inBound=.false.;ix_onLowBound=.false.;ix_onHighbound=.false.
 iy_inBound=.false.;iy_onLowBound=.false.;iy_onHighbound=.false.
 ir_inBound=.false.;ir_onTopBound=.false.;ir_onBotBound=.false.

 if(ixele.le.Tele_ixHighSub.and.ixele.ge.Tele_ixLowSub) then
    ix_inBound=.true.
 end if
 if(ixele.eq.Tele_ixLowSub) then
    ix_onLowBound=.true.
 end if
 if(ixele.eq.Tele_ixHighSub) then
    ix_onHighbound=.true.
 end if
 if(iyele.le.Tele_iyHighSub.and.iyele.ge.Tele_iyLowSub) then
     iy_inBound=.true.
 end if
 if(iyele.eq.Tele_iyLowSub) then
    iy_onLowBound=.true.
 end if
 if(iyele.eq.Tele_iyHighSub) then
    iy_onHighBound=.true.
 end if

 !if(isubregion.eq.1) then
 !   irele_true=irele
 !else if(isubregion.eq.2.or.isubregion.eq.3) then
 !   irele_true=irele+1 
 !else if(isubregion.eq.4.or.isubregion.eq.5) then
 !   irele_true=irele+2
 !end if
 
!The below lines are modified to match the new version SPECFEM3D. 
  if(isubregion.eq.1.or.isubregion.eq.2) then
     irele_true=irele     
  else if(isubregion.eq.3.or.isubregion.eq.4) then
!     irele_true=irele+1 
     irele_true=irele     
  else if(isubregion.eq.5) then
!     irele_true=irele+2
     irele_true=irele     
  end if
 


 if(irele_true.ge.Tele_irSub_Bot.and.irele_true.le.Tele_irSub_Top) then
    ir_inBound=.true.
 end if
 if(irele_true.eq.Tele_irSub_Top) then
    ir_onTopBound=.true.
 end if
 if(irele_true.eq.Tele_irSub_Bot) then
    ir_onBotBound=.true.
 end if

!debug
!if(myrank.eq.0) print *,'check2 ix_low',ixele,iyele,irele,Tele_ixHighSub,isubregion

  if(ix_onLowBound.and.ir_inBound.and.iy_inBound) then
     if(isubregion.eq.2.or.isubregion.eq.4) then
       if(ispec_superbrick.eq.17.or.ispec_superbrick.eq.18.or.ispec_superbrick.eq.21.or.&
          ispec_superbrick.eq.24.or.ispec_superbrick.eq.25.or.ispec_superbrick.eq.26.or. &
          ispec_superbrick.eq.29.or.ispec_superbrick.eq.32) then
            TeleEle_ixLow=TeleEle_ixLow+1

            TeleEle_xLow(TeleEle_ixLow)=ispec
            TeleEle_xLowReg(TeleEle_ixLow)=isubregion
!dubug
!            if(myrank.eq.0) print *,'check1 ix_low',ixele,iyele,irele,Tele_ixHighSub,isubregion
        end if
     else
        TeleEle_ixLow=TeleEle_ixLow+1
        TeleEle_xLow(TeleEle_ixLow)=ispec
        TeleEle_xLowReg(TeleEle_ixLow)=isubregion
!debug
!        if(myrank.eq.0) print *,'check ix_low',ixele,iyele,irele,Tele_ixHighSub,isubregion,&
!                                Tele_irSub_Top,Tele_irSub_Bot
     end if
     if(TeleEle_ixLow>TeleEle_nxLow)  then
      print   *,'ixele',ixele,iyele,'irele',irele,Tele_ixLowSub,&
                 Tele_ixHighSub,Tele_iyLowSub,Tele_iyHighSub,Tele_irSub_Top,Tele_irSub_Bot
      print *,'Error,The number of x-lower boundary &
              &elements for Teleseis-SEM outside the pridicted number',&
              TeleEle_ixLow,TeleEle_nxLow,iproc_xi,iproc_eta
     end if

  end if



  if(ix_onHighBound.and.ir_inBound.and.iy_inBound) then
     if(isubregion.eq.2.or.isubregion.eq.4) then
       if(ispec_superbrick.eq.1.or.ispec_superbrick.eq.2.or.ispec_superbrick.eq.5.or.&
          ispec_superbrick.eq.8.or.ispec_superbrick.eq.9.or.ispec_superbrick.eq.10.or. &
          ispec_superbrick.eq.13.or.ispec_superbrick.eq.16) then
            TeleEle_ixHigh=TeleEle_ixHigh +1
            if(TeleEle_ixHigh>TeleEle_nxHigh) stop 'Error,The number of x-higher boundary  &
            &elements for Teleseis-SEM outside the pridicted number'
            TeleEle_xHigh(TeleEle_ixHigh)=ispec
            TeleEle_xHighReg(TeleEle_ixHigh)=isubregion
        end if

     else
        TeleEle_ixHigh=TeleEle_ixHigh+1
        if(TeleEle_ixHigh>TeleEle_nxHigh)then
             print *,'ixele',ixele,iyele,'irele',irele,Tele_ixLowSub,Tele_ixHighSub,Tele_iyLowSub,Tele_iyHighSub,&
                      Tele_irSub_Top,Tele_irSub_Bot,TeleEle_ixHigh,TeleEle_nxHigh
        end if
        if(TeleEle_ixHigh>TeleEle_nxHigh) stop 'Error,The number of x-higher boundary  &
        &elements for Teleseis-SEM outside the pridicted number'

        TeleEle_xHigh(TeleEle_ixHigh)=ispec
        TeleEle_xHighReg(TeleEle_ixHigh)=isubregion
     end if
  end if


  if(iy_onLowBound.and.ir_inBound.and.ix_inBound) then
     if(isubregion.eq.2.or.isubregion.eq.4) then
        if(ispec_superbrick.eq.10.or.ispec_superbrick.eq.11.or.ispec_superbrick.eq.14.or.&
           ispec_superbrick.eq.16.or.ispec_superbrick.eq.26.or.ispec_superbrick.eq.27.or. &
           ispec_superbrick.eq.30.or.ispec_superbrick.eq.32) then
             TeleEle_iyLow=TeleEle_iyLow+1
!debug
!            print *,'iy-coupling',TeleEle_iyLow,TeleEle_nyLow,iyele,Tele_iyLowSub,myrank

             if(TeleEle_iyLow>TeleEle_nyLow) stop 'Error-superbrick,The number of y-lower boundary  &
             &elements for Teleseis-SEM outside the pridicted number'

             TeleEle_yLow(TeleEle_iyLow)=ispec
             TeleEle_yLowReg(TeleEle_iyLow)=isubregion
         end if
     else
        TeleEle_iyLow=TeleEle_iyLow+1
        if(TeleEle_iyLow>TeleEle_nyLow) stop 'Error,The number of y-lower boundary  &
        &elements for Teleseis-SEM outside the pridicted number'

        TeleEle_yLow(TeleEle_iyLow)=ispec
        TeleEle_yLowReg(TeleEle_iyLow)=isubregion
      end if
  end if


  if(iy_onHighBound.and.ir_inBound.and.ix_inBound) then
     if(isubregion.eq.2.or.isubregion.eq.4) then
        if(ispec_superbrick.eq.2.or.ispec_superbrick.eq.3.or.ispec_superbrick.eq.6.or. &
           ispec_superbrick.eq.8.or.ispec_superbrick.eq.18.or.ispec_superbrick.eq.19.or. &
           ispec_superbrick.eq.22.or.ispec_superbrick.eq.23) then
            TeleEle_iyHigh=TeleEle_iyHigh+1
            if(TeleEle_iyHigh>TeleEle_nyHigh) stop 'Error,The number of y-Higher boundary  &
            &elements for Teleseis-SEM outside the pridicted number'

            TeleEle_yHigh(TeleEle_iyHigh)=ispec
            TeleEle_yHighReg(TeleEle_iyHigh)=isubregion
         end if
     else
        TeleEle_iyHigh=TeleEle_iyHigh+1
        if(TeleEle_iyHigh>TeleEle_nyHigh) stop 'Error,The number of y-higher boundary  &
        &elements for Teleseis-SEM outside the pridicted number'

        TeleEle_yHigh(TeleEle_iyHigh)=ispec
        TeleEle_yHighReg(TeleEle_iyHigh)=isubregion
     end if
     if(TeleEle_iyHigh>TeleEle_nyHigh) then
         print *,'Error,The number of y-Higher boundary elements for &
                  &Teleseis-SEM outside the pridicted number'
         print *,'ixele',ixele,iyele,'irele',irele
     end if
  end if


  if(ir_onTopBound.and.ix_inBound.and.iy_inBound) then
      if(isubregion.eq.2.or.isubregion.eq.4) then
         if(ispec_superbrick.eq.8.or.ispec_superbrick.eq.16.or. &
            ispec_superbrick.eq.24.or.ispec_superbrick.eq.32) then
              TeleEle_irTop=TeleEle_irTop+1
!              if(TeleEle_ir>TeleEle_nr) stop 'Error,The number of r-bottom boundary  &
!              elements for Teleseis-SEM outside the pridicted number'

              TeleEle_rTop(TeleEle_irTop)=ispec
              TeleEle_rTopReg(TeleEle_irTop)=isubregion
         end if
      else
         TeleEle_irTop=TeleEle_irTop+1
!         if(TeleEle_ir>TeleEle_nr) stop 'Error,The number of r-bottom boundary  &
!         elements for Teleseis-SEM outside the pridicted number'

         TeleEle_rTop(TeleEle_irTop)=ispec
         TeleEle_rTopReg(TeleEle_irTop)=isubregion
      end if
      if(TeleEle_iyHigh>TeleEle_nyHigh) then
           !print *,'Error,The number of  r-bottom boundary elements for Teleseis-SEM outside the pridicted number'
           !print *,'ixele',ixele,iyele,'irele',irele
           !print *,'TeleEle_iyHigh',TeleEle_iyHigh,TeleEle_nyHigh
      end if
  end if
!  print *,'End getEle',myrank
  if(ir_onBotBound.and.ix_inBound.and.iy_inBound) then
      if(isubregion.eq.2.or.isubregion.eq.4) then
         if(ispec_superbrick.eq.8.or.ispec_superbrick.eq.16.or. &
            ispec_superbrick.eq.24.or.ispec_superbrick.eq.32) then
              TeleEle_irBot=TeleEle_irBot+1
!              if(TeleEle_ir>TeleEle_nr) stop 'Error,The number of r-bottom
!              boundary elements for Teleseis-SEM outside the pridicted number'

              TeleEle_rBot(TeleEle_irBot)=ispec
              TeleEle_rBotReg(TeleEle_irBot)=isubregion
         end if
      else
         TeleEle_irBot=TeleEle_irBot+1
!         if(TeleEle_ir>TeleEle_nr) stop 'Error,The number of r-bottom boundary&
!         elements for Teleseis-SEM outside the pridicted number'

         TeleEle_rBot(TeleEle_irBot)=ispec
         TeleEle_rBotReg(TeleEle_irBot)=isubregion
      end if
      if(TeleEle_iyHigh>TeleEle_nyHigh) then
           !print *,'Error,The number of  r-bottom boundary elements for Teleseis-SEM outside the pridicted number'
           !print *,'ixele',ixele,iyele,'irele',irele
           !print *,'TeleEle_iyHigh',TeleEle_iyHigh,TeleEle_nyHigh
      end if
  end if

end subroutine getEle_TeleBound

             
