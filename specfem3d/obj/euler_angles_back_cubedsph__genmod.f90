        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:59 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EULER_ANGLES_BACK_CUBEDSPH__genmod
          INTERFACE 
            SUBROUTINE EULER_ANGLES_BACK_CUBEDSPH(                      &
     &ROTATION_MATRIX_BACK_CUBEDSPH,CENTER_LONGITUDE_IN_DEGREES,        &
     &CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)
              REAL(KIND=8) :: ROTATION_MATRIX_BACK_CUBEDSPH(3,3)
              REAL(KIND=8) :: CENTER_LONGITUDE_IN_DEGREES
              REAL(KIND=8) :: CENTER_LATITUDE_IN_DEGREES
              REAL(KIND=8) :: GAMMA_ROTATION_AZIMUTH
            END SUBROUTINE EULER_ANGLES_BACK_CUBEDSPH
          END INTERFACE 
        END MODULE EULER_ANGLES_BACK_CUBEDSPH__genmod
