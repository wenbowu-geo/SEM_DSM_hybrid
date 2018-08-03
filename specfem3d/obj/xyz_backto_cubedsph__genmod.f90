        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:08 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE XYZ_BACKTO_CUBEDSPH__genmod
          INTERFACE 
            SUBROUTINE XYZ_BACKTO_CUBEDSPH(X,Y,Z,X_CUBEDSPH,Y_CUBEDSPH, &
     &Z_CUBEDSPH,ANGULAR_WIDTH_XI_IN_DEGREES,                           &
     &ANGULAR_WIDTH_ETA_IN_DEGREES,ROTATION_MATRIX_BACK_CUBEDSPH)
              REAL(KIND=8), INTENT(IN) :: X
              REAL(KIND=8), INTENT(IN) :: Y
              REAL(KIND=8), INTENT(IN) :: Z
              REAL(KIND=8), INTENT(OUT) :: X_CUBEDSPH
              REAL(KIND=8), INTENT(OUT) :: Y_CUBEDSPH
              REAL(KIND=8), INTENT(OUT) :: Z_CUBEDSPH
              REAL(KIND=8) :: ANGULAR_WIDTH_XI_IN_DEGREES
              REAL(KIND=8) :: ANGULAR_WIDTH_ETA_IN_DEGREES
              REAL(KIND=8) :: ROTATION_MATRIX_BACK_CUBEDSPH(3,3)
            END SUBROUTINE XYZ_BACKTO_CUBEDSPH
          END INTERFACE 
        END MODULE XYZ_BACKTO_CUBEDSPH__genmod
