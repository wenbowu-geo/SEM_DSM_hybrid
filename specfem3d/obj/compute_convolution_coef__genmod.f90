        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_CONVOLUTION_COEF__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_CONVOLUTION_COEF(BB,DELTAT,COEF0,COEF1,  &
     &COEF2)
              REAL(KIND=8), INTENT(IN) :: BB
              REAL(KIND=8), INTENT(IN) :: DELTAT
              REAL(KIND=8), INTENT(OUT) :: COEF0
              REAL(KIND=8), INTENT(OUT) :: COEF1
              REAL(KIND=8), INTENT(OUT) :: COEF2
            END SUBROUTINE COMPUTE_CONVOLUTION_COEF
          END INTERFACE 
        END MODULE COMPUTE_CONVOLUTION_COEF__genmod
