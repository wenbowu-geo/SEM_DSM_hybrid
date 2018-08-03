!**********************************allocate the struct size*************************************************
!subroutine Receiver_default(ReceiverInfo)
subroutine Receiver_default()
   use SEMtoTele_par
   implicit none
   integer ::i
!Type(Receiver) ::ReceiverInfo
   !   ReceiverInfo%lamdalReveiver=-9999.0
   !   ReceiverInfo%muReceiver   =-9999.0
   !   ReceiverInfo%VpReceiver    =-9999.0
   !   ReceiverInfo%RhoReceiver   =-9999.0
  do i=1,Nrecv
   ReceiverInfo(i)%VpReceiver    =VpReceiver
   ReceiverInfo(i)%RhoReceiver   =RhoReceiver
   ReceiverInfo(i)%muReceiver    =VsReceiver**2*RhoReceiver/2.0
   ReceiverInfo(i)%lamdalReveiver=(VpReceiver**2-VsReceiver**2)*RhoReceiver
   !ReceiverInfo(i)%eRadi(:)      =0.0
   !ReceiverInfo(i)%eNorth(:)     =0.0
   !ReceiverInfo(i)%eEast(:)      =0.0
   ReceiverInfo(i)%lat           =-9999.0
   ReceiverInfo(i)%lon           =-9999.0
   ReceiverInfo(i)%Theta         =-9999.0
   ReceiverInfo(i)%Phi           =-9999.0
   ReceiverInfo(i)%Radi          =-1.e10
   ReceiverInfo(i)%x             =-1.e10
   ReceiverInfo(i)%y             =-1.e10
   ReceiverInfo(i)%z             =-1.e10
  end do
end subroutine Receiver_default


!subroutine SeismoRec_default(Seismogram,DT,NSTEP,myrank)
!use struct_defined
!use constant_SEMtoTele
!implicit none
!Type(SeismoRec) ::Seismogram
!double precision::DT
!integer         ::NSTEP
!integer         ::myrank
!Seismogram%dt=DT
!Seismogram%Nt=NTMAX
!print *,'Seismogram1=',Seismogram%dt,DT,Seismogram%Nt,NTMAX,myrank
!allocate(Seismogram%Seismo(3,NTMAX))
!Seismogram%Seismo(:,:)=0.0
!print *,'Seismogram=',Seismogram%dt,DT,Seismogram%Nt,NTMAX,myrank
!end subroutine SeismoRec_default

subroutine SeismoRec_default()
use SEMtoTele_par
use specfem_par
implicit none
integer ::i
 do i=1,NRecv
   Seismogram(i)%dt=DT
   Seismogram(i)%nt=NTMAX
!   allocate(Seismogram(i)%Seismo(3,NTMAX))
   Seismogram(i)%Seismo(:,:)=0.0
 end do
end subroutine SeismoRec_default


!subroutine ViaRepresent_default(RepInfo,NGLLX,NGLLY,NGLLZ,Nele_Bound,NRecv)
!subroutine ViaRepresent_default()
!    use SEMtoTele_par
!   Type(ViaRepresent)::RepInfo
!   integer           ::NGLLX,NGLLY,NGLLZ,Nele_Bound,NRecv
!   allocate(RepInfo%ispec_bound(Nele_Bound))
!   allocate(RepInfo%iregion(Nele_Bound))
!   allocate(RepInfo%type_Bound(Nele_Bound))
!   allocate(RepInfo%x(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%y(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%z(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%Theta(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%Phi(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%Radi(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%velocity(3,NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%sigmaxx(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%sigmaxy(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%sigmaxz(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%sigmayy(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%sigmayx(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%sigmayz(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%sigmazx(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%sigmazy(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%sigmazz(NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%Travel_time(NGLLX,NGLLY,NGLLZ,Nele_Bound,NRecv))
!   allocate(RepInfo%RayVectorReceiver(3,NGLLX,NGLLY,NGLLZ,Nele_Bound,NRecv))
!   allocate(RepInfo%RayVectorBound(3,NGLLX,NGLLY,NGLLZ,Nele_Bound,NRecv))
!   allocate(RepInfo%VectorNormal(3,NGLLX,NGLLY,NGLLZ,Nele_Bound))
!   allocate(RepInfo%Gcarc(NGLLX,NGLLY,NGLLZ,Nele_Bound,NRecv))
!   allocate(RepInfo%Geofactor(NGLLX,NGLLY,NGLLZ,Nele_Bound,NRecv))
!end subroutine
