        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:32 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ATTENUATION_CONSTANTS__genmod
          INTERFACE 
            SUBROUTINE GET_ATTENUATION_CONSTANTS(MIN_RESOLVED_PERIOD,   &
     &TAU_SIGMA,F_C_SOURCE,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD&
     &)
              REAL(KIND=8) :: MIN_RESOLVED_PERIOD
              REAL(KIND=8) :: TAU_SIGMA(3)
              REAL(KIND=8) :: F_C_SOURCE
              REAL(KIND=8) :: MIN_ATTENUATION_PERIOD
              REAL(KIND=8) :: MAX_ATTENUATION_PERIOD
            END SUBROUTINE GET_ATTENUATION_CONSTANTS
          END INTERFACE 
        END MODULE GET_ATTENUATION_CONSTANTS__genmod
