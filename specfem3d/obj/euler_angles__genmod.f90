        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:03 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EULER_ANGLES__genmod
          INTERFACE 
            SUBROUTINE EULER_ANGLES(ROTATION_MATRIX,                    &
     &CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,           &
     &GAMMA_ROTATION_AZIMUTH)
              REAL(KIND=8) :: ROTATION_MATRIX(3,3)
              REAL(KIND=8) :: CENTER_LONGITUDE_IN_DEGREES
              REAL(KIND=8) :: CENTER_LATITUDE_IN_DEGREES
              REAL(KIND=8) :: GAMMA_ROTATION_AZIMUTH
            END SUBROUTINE EULER_ANGLES
          END INTERFACE 
        END MODULE EULER_ANGLES__genmod
