        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:19 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LX_PARAMETER_COMPUTATION__genmod
          INTERFACE 
            SUBROUTINE LX_PARAMETER_COMPUTATION(DELTAT,KAPPA_X,D_X,     &
     &ALPHA_X,CPML_REGION_LOCAL,A_0,A_1,COEF0_1,COEF1_1,COEF2_1)
              REAL(KIND=8), INTENT(IN) :: DELTAT
              REAL(KIND=8), INTENT(IN) :: KAPPA_X
              REAL(KIND=8), INTENT(IN) :: D_X
              REAL(KIND=8), INTENT(IN) :: ALPHA_X
              INTEGER(KIND=4), INTENT(IN) :: CPML_REGION_LOCAL
              REAL(KIND=8), INTENT(OUT) :: A_0
              REAL(KIND=8), INTENT(OUT) :: A_1
              REAL(KIND=8), INTENT(OUT) :: COEF0_1
              REAL(KIND=8), INTENT(OUT) :: COEF1_1
              REAL(KIND=8), INTENT(OUT) :: COEF2_1
            END SUBROUTINE LX_PARAMETER_COMPUTATION
          END INTERFACE 
        END MODULE LX_PARAMETER_COMPUTATION__genmod
