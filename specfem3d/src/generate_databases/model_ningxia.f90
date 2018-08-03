 subroutine model_ningxia(x_eval,y_eval,z_eval,r_middle, &
                             rho_final,vp_final,vs_final,qkappa_atten,qmu_atten,imaterial_id)

  use tomography

  implicit none

  double precision, intent(in) :: x_eval,y_eval,z_eval,r_middle
  real(kind=CUSTOM_REAL), intent(out) :: vp_final,vs_final,rho_final
  real(kind=CUSTOM_REAL), intent(out) :: qkappa_atten,qmu_atten
  integer, intent(in) :: imaterial_id

  integer ::myrank

  ! local parameters

  double precision ::rho_layer(4,3),vp_layer(4,3),vs_layer(4,3)
  double precision ::r_layer(3)
  double precision ::rearth,Radi,rend,temp
  integer ::ilayer

  qmu_atten = 80.0
  qkappa_atten=9999.

  if(imaterial_id.eq.1) then
     vp_final=6500.0
     vs_final=3700.0
     rho_final=2900.0
  else if(imaterial_id.eq.2) then
     vp_final=8000.0
     vs_final=4700.0
     rho_final=2380.0
  else 
     stop 'Error tomography_ningxia'
  end if



  call get_perturbation(x_eval,y_eval,z_eval, &
                             rho_final,vp_final,vs_final)
  !print
  !*,'structure',6371000.0-Radi,vp_final,vs_final,rho_final,x_eval,y_eval,z_eval
!  if(vs_final.lt.1.0) stop 'vs error'


  end subroutine model_ningxia

