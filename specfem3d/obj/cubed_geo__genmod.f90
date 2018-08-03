        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:43 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CUBED_GEO__genmod
          INTERFACE 
            SUBROUTINE CUBED_GEO(XELM,YELM,ZELM,                        &
     &ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,         &
     &ROTATION_MATRIX)
              REAL(KIND=8) :: XELM(8)
              REAL(KIND=8) :: YELM(8)
              REAL(KIND=8) :: ZELM(8)
              REAL(KIND=8) :: ANGULAR_WIDTH_XI_IN_DEGREES
              REAL(KIND=8) :: ANGULAR_WIDTH_ETA_IN_DEGREES
              REAL(KIND=8) :: ROTATION_MATRIX(3,3)
            END SUBROUTINE CUBED_GEO
          END INTERFACE 
        END MODULE CUBED_GEO__genmod
