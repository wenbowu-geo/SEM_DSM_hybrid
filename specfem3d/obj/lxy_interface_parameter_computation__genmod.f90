        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:19 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LXY_INTERFACE_PARAMETER_COMPUTATION__genmod
          INTERFACE 
            SUBROUTINE LXY_INTERFACE_PARAMETER_COMPUTATION(DELTAT,      &
     &KAPPA_X,D_X,ALPHA_X,KAPPA_Y,D_Y,ALPHA_Y,CPML_REGION_LOCAL,        &
     &INDEX_IJK,A_0,A_1,A_2,COEF0_1,COEF1_1,COEF2_1,COEF0_2,COEF1_2,    &
     &COEF2_2)
              REAL(KIND=8), INTENT(IN) :: DELTAT
              REAL(KIND=8) :: KAPPA_X
              REAL(KIND=8) :: D_X
              REAL(KIND=8) :: ALPHA_X
              REAL(KIND=8) :: KAPPA_Y
              REAL(KIND=8) :: D_Y
              REAL(KIND=8) :: ALPHA_Y
              INTEGER(KIND=4), INTENT(IN) :: CPML_REGION_LOCAL
              INTEGER(KIND=4), INTENT(IN) :: INDEX_IJK
              REAL(KIND=8), INTENT(OUT) :: A_0
              REAL(KIND=8), INTENT(OUT) :: A_1
              REAL(KIND=8), INTENT(OUT) :: A_2
              REAL(KIND=8), INTENT(OUT) :: COEF0_1
              REAL(KIND=8), INTENT(OUT) :: COEF1_1
              REAL(KIND=8), INTENT(OUT) :: COEF2_1
              REAL(KIND=8), INTENT(OUT) :: COEF0_2
              REAL(KIND=8), INTENT(OUT) :: COEF1_2
              REAL(KIND=8), INTENT(OUT) :: COEF2_2
            END SUBROUTINE LXY_INTERFACE_PARAMETER_COMPUTATION
          END INTERFACE 
        END MODULE LXY_INTERFACE_PARAMETER_COMPUTATION__genmod
