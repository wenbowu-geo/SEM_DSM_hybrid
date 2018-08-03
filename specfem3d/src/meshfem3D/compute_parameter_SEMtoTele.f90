subroutine compute_parameter_SEMtoTele(Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,&
                                 Tele_irTop,Tele_irBot,TeleEle_nxLow, &
                                 TeleEle_nxHigh,TeleEle_nyLow,TeleEle_nyHigh,TeleEle_nrTop,&
                                 TeleEle_nrBot,iproc_xi,iproc_eta, &
                                 NEX_XI,NEX_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER, &
                                 USE_REGULAR_MESH,NDOUBLINGS,ner_doublings,DEBUG_COUPLING,myrank)

  implicit none
  logical USE_REGULAR_MESH
  integer myrank
  logical ::DEBUG_COUPLING
  integer NDOUBLINGS
  integer:: Tele_ixLow,Tele_ixHigh,Tele_iyLow,Tele_iyHigh,Tele_irTop,Tele_irBot
  integer::TeleEle_nxLow,TeleEle_nxHigh,TeleEle_nyLow,TeleEle_nyHigh,TeleEle_nrTop,TeleEle_nrBot
  integer::iproc_xi,iproc_eta
  integer::NEX_XI,NEX_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER
  integer, dimension(2) :: ner_doublings
!  integer nblayers
!  integer ner_layer(nblayers)
  
  integer::ixMin,ixMax,iyMin,iyMax
  logical::ixLow_inside,ixHigh_inside,Ele_insidexLowHigh,iyLow_inside,iyHigh_inside,Ele_insideyLowHigh
  integer::dix,diy,dix1,diy1,dir1,dix2,diy2,dir2,dix3,diy3,dir3,dix4,diy4,dir4,dix5,diy5,dir5

 integer ::int_tmp

!The below three lines are used to avoid possible error reports during compling
!the code. NER and myrank are not used now, that occationally causes error 
!reports. But they might be usful in the future, so we keep them here.
 int_tmp=NEX_XI
 int_tmp=NEX_ETA
 int_tmp=NER


  TeleEle_nxLow=0;TeleEle_nxHigh=0
  TeleEle_nyLow=0;TeleEle_nyHigh=0
  TeleEle_nrTop=0;TeleEle_nrBot =0
  ixMin=iproc_xi*NEX_PER_PROC_XI+1
  ixMax=(iproc_xi+1)*NEX_PER_PROC_XI
  iyMin=iproc_eta*NEX_PER_PROC_ETA+1
  iyMax=(iproc_eta+1)*NEX_PER_PROC_ETA
  ixLow_inside=.false.;ixHigh_inside=.false.;Ele_insidexLowHigh=.false.
  iyLow_inside=.false.;iyHigh_inside=.false.;Ele_insideyLowHigh=.false.

  if(ixMin.le.Tele_ixLow.and.ixMax.ge.Tele_ixLow) then
      ixLow_inside=.true.
  endif

  if(ixMin.le.Tele_ixHigh.and.ixMax.ge.Tele_ixHigh) then
      ixHigh_inside=.true.
  end if
  if(ixMin.ge.Tele_ixLow.and.ixMax.le.Tele_ixHigh) then
      Ele_insidexLowHigh=.true.
  end if
  if(iyMin.le.Tele_iyLow.and.iyMax.ge.Tele_iyLow) then
       iyLow_inside=.true.
  end if
  if(iyMin.le.Tele_iyHigh.and.iyMax.ge.Tele_iyHigh) then
       iyHigh_inside=.true.
  end if
  if(iyMin.ge.Tele_iyLow.and.iyMax.le.Tele_iyHigh) then
       Ele_insideyLowHigh=.true.
  end if


  dix=0
  if(ixLow_inside.and.ixHigh_inside) then
      dix=Tele_ixHigh-Tele_ixLow+1
  else if(ixLow_inside.and.(.not.ixHigh_inside)) then
      dix=ixMax-Tele_ixLow+1
  else if(.not.ixLow_inside.and.ixHigh_inside) then
      dix=Tele_ixHigh-ixMin+1
  else if(Ele_insidexLowHigh) then
      dix=NEX_PER_PROC_XI
  end if

  diy=0
  if(iyLow_inside.and.iyHigh_inside) then
      diy=Tele_iyHigh-Tele_iyLow+1
  else if(iyLow_inside.and.(.not.iyHigh_inside)) then
      diy=iyMax-Tele_iyLow+1
  else if(.not.iyLow_inside.and.iyHigh_inside) then
      diy=Tele_iyHigh-iyMin+1
  else if(Ele_insideyLowHigh) then
      diy=NEX_PER_PROC_ETA
  end if


  if(USE_REGULAR_MESH) then
   dix1=dix;diy1=diy;dir1=Tele_irTop-Tele_irBot+1
   dix2=0;diy2=0;dir2=0
   dix3=0;diy3=0;dir3=0
   dix4=0;diy4=0;dir4=0
   dix5=0;diy5=0;dir5=0
   TeleEle_nrTop=dix1*diy1
   TeleEle_nrBot=dix1*diy1
  else if(NDOUBLINGS.eq.1) then
   dix1=dix/2;diy1=diy/2
   dix2=dix/4;diy2=diy/4
   dix3=dix;diy3=diy
   dix4=0;diy4=0
   dix5=0;diy5=0
   if(ner_doublings(1)-2.gt.Tele_irBot.and.ner_doublings(1)-2.lt.Tele_irTop) then
       dir1=ner_doublings(1)-2-Tele_irBot+1;dir2=2;dir3=Tele_irTop-ner_doublings(1);dir4=0;dir5=0
       TeleEle_nrBot=dix1*diy1
       TeleEle_nrTop=dix3*diy3
   else 
       stop 'Tele_irTop is below doubling layer or Tele_irBot is above it'
   end if
 !  else if(ner_doublings(1).eq.Tele_ir.or.ner_doublings(1)-1.eq.Tele_ir) then
 !      dir1=0;dir2=2;dir3=NER-ner_doublings(1);dir4=0;dir5=0
 !      TeleEle_nrTop=4*dix2*diy2
 !  else 
 !    dir1=0;dir2=0;dir3=NER-Tele_ir+1;dir4=0;dir5=0
 !      TeleEle_nr=dix3*diy3
 !  end if
!   print *,'TeleEle_nr=',TeleEle_nr,dix,diy
   
  else if(NDOUBLINGS.eq.2) then
   dix1=dix/4;diy1=diy/4
   dix2=dix/8;diy2=diy/8
   dix3=dix/2;diy3=diy/2
   dix4=dix/4;diy4=diy/4
   dix5=dix;diy5=diy
   if(ner_doublings(2)-2.gt.Tele_irBot.and.ner_doublings(1).lt.Tele_irTop) then
       dir1=ner_doublings(2)-2-Tele_irBot+1;dir2=2;dir3=(ner_doublings(1) - 2) - ner_doublings(2)
       dir4=2;dir5=Tele_irTop-ner_doublings(1)
       TeleEle_nrBot=dix1*diy1
       TeleEle_nrTop=dix5*diy5
   else 
       stop 'Tele_irTop is below doubling layer or Tele_irBot is above it'
!   else if (ner_doublings(2).eq.Tele_ir.or.ner_doublings(2)-1.eq.Tele_ir) then
!       dir1=0;dir2=2;dir3=(ner_doublings(1) - 2) - ner_doublings(2)
!       dir4=2;dir5=NER-ner_doublings(1)
!       TeleEle_nr=4*dix2*diy2
!   else if(ner_doublings(2).lt.Tele_ir.and.ner_doublings(1)-1.gt.Tele_ir) then
!       dir1=0;dir2=0;dir3=ner_doublings(1) - 2 -Tele_ir+1
!       dir4=2;dir5=NER-ner_doublings(1)
!       TeleEle_nr=dix3*diy3
!   else if(ner_doublings(1).eq.Tele_ir.or.ner_doublings(1)-1.eq.Tele_ir) then
!       dir1=0;dir2=0;dir3=0;dir4=2;dir5=NER-ner_doublings(1)
!       TeleEle_nr=4*dix4*diy4
!   else
!       dir1=0;dir2=0;dir3=0;dir4=0;dir5=NER-Tele_ir+1
!       TeleEle_nr=dix5*diy5
   end if
!   print *,'TeleEle_nr1=',TeleEle_nr,dix,diy
  end if



  if(ixLow_inside) then
     TeleEle_nxLow=diy1*dir1+4*diy2*dir2+diy3*dir3+4*diy4*dir4+diy5*dir5
     if(DEBUG_COUPLING) then
         print  *,'TeleEle_nxLow1',TeleEle_nxLow,diy1*dir1,4*diy2*dir2,&
                            diy3*dir3,diy3,dix3,4*diy4*dir4,diy5*dir5,Ele_insidexLowHigh,&
                        ixMin,Tele_ixLow,ixMax,Tele_ixHigh,myrank,dir1,diy1,dir2,diy2,dir3,diy3
     end if
  end if 
  if(ixHigh_inside) then
     TeleEle_nxHigh=diy1*dir1+4*diy2*dir2+diy3*dir3+4*diy4*dir4+diy5*dir5
     if(DEBUG_COUPLING) then
         print  *,'TeleEle_nxHigh1',TeleEle_nxLow,diy1*dir1,4*diy2*dir2,diy3*dir3,4*diy4*dir4,diy5*dir5,myrank
     end if
  end if
  if(iyLow_inside) then
     TeleEle_nyLow= dix1*dir1+4*dix2*dir2+dix3*dir3+4*dix4*dir4+dix5*dir5
     if(DEBUG_COUPLING) then
         print     *,'TeleEle_nyLow1',dix1*dir1,4*dix2*dir2,dix3*dir3,4*dix4*dir4,dix5*dir5,myrank
     end if
  end if
  if(iyHigh_inside) then
     TeleEle_nyHigh=dix1*dir1+4*dix2*dir2+dix3*dir3+4*dix4*dir4+dix5*dir5
     if(DEBUG_COUPLING) then
        print     *,'TeleEle_nyHigh1',dix1*dir1,4*dix2*dir2,dix3*dir3,4*dix4*dir4,dix5*dir5,myrank
     end if
  end if
  if(DEBUG_COUPLING) then
        print *,'Nele_x_y_z',myrank,TeleEle_nxLow,TeleEle_nxHigh,TeleEle_nyLow,TeleEle_nyHigh,&
                TeleEle_nrBot,TeleEle_nrTop
  end if
end subroutine compute_parameter_SEMtoTele
