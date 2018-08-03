        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:56 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MODEL_TOMOGRAPHY1__genmod
          INTERFACE 
            SUBROUTINE MODEL_TOMOGRAPHY1(X_EVAL,Y_EVAL,Z_EVAL,R_MIDDLE, &
     &RHO_FINAL,VP_FINAL,VS_FINAL,MYRANK)
              REAL(KIND=8), INTENT(IN) :: X_EVAL
              REAL(KIND=8), INTENT(IN) :: Y_EVAL
              REAL(KIND=8), INTENT(IN) :: Z_EVAL
              REAL(KIND=8), INTENT(IN) :: R_MIDDLE
              REAL(KIND=8), INTENT(OUT) :: RHO_FINAL
              REAL(KIND=8), INTENT(OUT) :: VP_FINAL
              REAL(KIND=8), INTENT(OUT) :: VS_FINAL
              INTEGER(KIND=4) :: MYRANK
            END SUBROUTINE MODEL_TOMOGRAPHY1
          END INTERFACE 
        END MODULE MODEL_TOMOGRAPHY1__genmod
