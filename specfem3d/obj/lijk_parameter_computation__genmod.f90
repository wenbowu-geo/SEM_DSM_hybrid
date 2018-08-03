        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:19 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LIJK_PARAMETER_COMPUTATION__genmod
          INTERFACE 
            SUBROUTINE LIJK_PARAMETER_COMPUTATION(DELTAT,KAPPA_X,D_X,   &
     &ALPHA_X,KAPPA_Y,D_Y,ALPHA_Y,KAPPA_Z,D_Z,ALPHA_Z,CPML_REGION_LOCAL,&
     &INDEX_IJK,A_0,A_1,A_2,A_3,COEF0_1,COEF1_1,COEF2_1,COEF0_2,COEF1_2,&
     &COEF2_2,COEF0_3,COEF1_3,COEF2_3)
              REAL(KIND=8), INTENT(IN) :: DELTAT
              REAL(KIND=8), INTENT(IN) :: KAPPA_X
              REAL(KIND=8), INTENT(IN) :: D_X
              REAL(KIND=8), INTENT(IN) :: ALPHA_X
              REAL(KIND=8), INTENT(IN) :: KAPPA_Y
              REAL(KIND=8), INTENT(IN) :: D_Y
              REAL(KIND=8), INTENT(IN) :: ALPHA_Y
              REAL(KIND=8), INTENT(IN) :: KAPPA_Z
              REAL(KIND=8), INTENT(IN) :: D_Z
              REAL(KIND=8), INTENT(IN) :: ALPHA_Z
              INTEGER(KIND=4), INTENT(IN) :: CPML_REGION_LOCAL
              INTEGER(KIND=4), INTENT(IN) :: INDEX_IJK
              REAL(KIND=8), INTENT(OUT) :: A_0
              REAL(KIND=8), INTENT(OUT) :: A_1
              REAL(KIND=8), INTENT(OUT) :: A_2
              REAL(KIND=8), INTENT(OUT) :: A_3
              REAL(KIND=8), INTENT(OUT) :: COEF0_1
              REAL(KIND=8), INTENT(OUT) :: COEF1_1
              REAL(KIND=8), INTENT(OUT) :: COEF2_1
              REAL(KIND=8), INTENT(OUT) :: COEF0_2
              REAL(KIND=8), INTENT(OUT) :: COEF1_2
              REAL(KIND=8), INTENT(OUT) :: COEF2_2
              REAL(KIND=8), INTENT(OUT) :: COEF0_3
              REAL(KIND=8), INTENT(OUT) :: COEF1_3
              REAL(KIND=8), INTENT(OUT) :: COEF2_3
            END SUBROUTINE LIJK_PARAMETER_COMPUTATION
          END INTERFACE 
        END MODULE LIJK_PARAMETER_COMPUTATION__genmod
