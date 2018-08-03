        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE L_PARAMETER_COMPUTATION__genmod
          INTERFACE 
            SUBROUTINE L_PARAMETER_COMPUTATION(DELTAT,KAPPA_X,D_X,      &
     &ALPHA_X,KAPPA_Y,D_Y,ALPHA_Y,KAPPA_Z,D_Z,ALPHA_Z,CPML_REGION_LOCAL,&
     &A_0,A_1,A_2,A_3,A_4,A_5,COEF0_X,COEF1_X,COEF2_X,COEF0_Y,COEF1_Y,  &
     &COEF2_Y,COEF0_Z,COEF1_Z,COEF2_Z)
              REAL(KIND=8), INTENT(IN) :: DELTAT
              REAL(KIND=8) :: KAPPA_X
              REAL(KIND=8) :: D_X
              REAL(KIND=8) :: ALPHA_X
              REAL(KIND=8) :: KAPPA_Y
              REAL(KIND=8) :: D_Y
              REAL(KIND=8) :: ALPHA_Y
              REAL(KIND=8) :: KAPPA_Z
              REAL(KIND=8) :: D_Z
              REAL(KIND=8) :: ALPHA_Z
              INTEGER(KIND=4), INTENT(IN) :: CPML_REGION_LOCAL
              REAL(KIND=8), INTENT(OUT) :: A_0
              REAL(KIND=8), INTENT(OUT) :: A_1
              REAL(KIND=8), INTENT(OUT) :: A_2
              REAL(KIND=8), INTENT(OUT) :: A_3
              REAL(KIND=8), INTENT(OUT) :: A_4
              REAL(KIND=8), INTENT(OUT) :: A_5
              REAL(KIND=8), INTENT(OUT) :: COEF0_X
              REAL(KIND=8), INTENT(OUT) :: COEF1_X
              REAL(KIND=8), INTENT(OUT) :: COEF2_X
              REAL(KIND=8), INTENT(OUT) :: COEF0_Y
              REAL(KIND=8), INTENT(OUT) :: COEF1_Y
              REAL(KIND=8), INTENT(OUT) :: COEF2_Y
              REAL(KIND=8), INTENT(OUT) :: COEF0_Z
              REAL(KIND=8), INTENT(OUT) :: COEF1_Z
              REAL(KIND=8), INTENT(OUT) :: COEF2_Z
            END SUBROUTINE L_PARAMETER_COMPUTATION
          END INTERFACE 
        END MODULE L_PARAMETER_COMPUTATION__genmod
