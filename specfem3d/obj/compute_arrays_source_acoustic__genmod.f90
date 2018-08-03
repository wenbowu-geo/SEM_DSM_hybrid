        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:42 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_ARRAYS_SOURCE_ACOUSTIC__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_ARRAYS_SOURCE_ACOUSTIC(SOURCEARRAY,HXIS, &
     &HETAS,HGAMMAS,FACTOR_SOURCE)
              REAL(KIND=8) :: SOURCEARRAY(3,5,5,5)
              REAL(KIND=8) :: HXIS(5)
              REAL(KIND=8) :: HETAS(5)
              REAL(KIND=8) :: HGAMMAS(5)
              REAL(KIND=8) :: FACTOR_SOURCE
            END SUBROUTINE COMPUTE_ARRAYS_SOURCE_ACOUSTIC
          END INTERFACE 
        END MODULE COMPUTE_ARRAYS_SOURCE_ACOUSTIC__genmod
