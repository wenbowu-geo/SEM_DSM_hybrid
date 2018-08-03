        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:56 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MODEL_NINGXIA__genmod
          INTERFACE 
            SUBROUTINE MODEL_NINGXIA(X_EVAL,Y_EVAL,Z_EVAL,R_MIDDLE,     &
     &RHO_FINAL,VP_FINAL,VS_FINAL,QKAPPA_ATTEN,QMU_ATTEN,IMATERIAL_ID)
              REAL(KIND=8), INTENT(IN) :: X_EVAL
              REAL(KIND=8), INTENT(IN) :: Y_EVAL
              REAL(KIND=8), INTENT(IN) :: Z_EVAL
              REAL(KIND=8), INTENT(IN) :: R_MIDDLE
              REAL(KIND=8), INTENT(OUT) :: RHO_FINAL
              REAL(KIND=8), INTENT(OUT) :: VP_FINAL
              REAL(KIND=8), INTENT(OUT) :: VS_FINAL
              REAL(KIND=8), INTENT(OUT) :: QKAPPA_ATTEN
              REAL(KIND=8), INTENT(OUT) :: QMU_ATTEN
              INTEGER(KIND=4), INTENT(IN) :: IMATERIAL_ID
            END SUBROUTINE MODEL_NINGXIA
          END INTERFACE 
        END MODULE MODEL_NINGXIA__genmod
