        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:56 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_PERTURBATION__genmod
          INTERFACE 
            SUBROUTINE GET_PERTURBATION(X_EVAL,Y_EVAL,Z_EVAL,RHO_FINAL, &
     &VP_FINAL,VS_FINAL)
              REAL(KIND=8), INTENT(IN) :: X_EVAL
              REAL(KIND=8), INTENT(IN) :: Y_EVAL
              REAL(KIND=8), INTENT(IN) :: Z_EVAL
              REAL(KIND=8), INTENT(OUT) :: RHO_FINAL
              REAL(KIND=8), INTENT(OUT) :: VP_FINAL
              REAL(KIND=8), INTENT(OUT) :: VS_FINAL
            END SUBROUTINE GET_PERTURBATION
          END INTERFACE 
        END MODULE GET_PERTURBATION__genmod
