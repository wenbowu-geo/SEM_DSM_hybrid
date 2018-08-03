        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:56 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MODEL_TOMOGRAPHY__genmod
          INTERFACE 
            SUBROUTINE MODEL_TOMOGRAPHY(FLAG_MEDIA,X_EVAL,Y_EVAL,Z_EVAL,&
     &XEVAL_CUBEDSPH,YEVAL_CUBEDSPH,ZEVAL_CUBEDSPH,R_MIDDLE,RHO_FINAL,  &
     &VP_FINAL,VS_FINAL,QKAPPA_ATTEN,QMU_ATTEN,MYRANK)
              INTEGER(KIND=4), INTENT(IN) :: FLAG_MEDIA
              REAL(KIND=8), INTENT(IN) :: X_EVAL
              REAL(KIND=8), INTENT(IN) :: Y_EVAL
              REAL(KIND=8), INTENT(IN) :: Z_EVAL
              REAL(KIND=8), INTENT(IN) :: XEVAL_CUBEDSPH
              REAL(KIND=8), INTENT(IN) :: YEVAL_CUBEDSPH
              REAL(KIND=8), INTENT(IN) :: ZEVAL_CUBEDSPH
              REAL(KIND=8), INTENT(IN) :: R_MIDDLE
              REAL(KIND=8), INTENT(OUT) :: RHO_FINAL
              REAL(KIND=8), INTENT(OUT) :: VP_FINAL
              REAL(KIND=8), INTENT(OUT) :: VS_FINAL
              REAL(KIND=8), INTENT(OUT) :: QKAPPA_ATTEN
              REAL(KIND=8), INTENT(OUT) :: QMU_ATTEN
              INTEGER(KIND=4), INTENT(IN) :: MYRANK
            END SUBROUTINE MODEL_TOMOGRAPHY
          END INTERFACE 
        END MODULE MODEL_TOMOGRAPHY__genmod
