        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:32 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ATTENUATION_PERIODS__genmod
          INTERFACE 
            SUBROUTINE GET_ATTENUATION_PERIODS(MIN_RESOLVED_PERIOD,     &
     &MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)
              REAL(KIND=8), INTENT(IN) :: MIN_RESOLVED_PERIOD
              REAL(KIND=8), INTENT(OUT) :: MIN_ATTENUATION_PERIOD
              REAL(KIND=8), INTENT(OUT) :: MAX_ATTENUATION_PERIOD
            END SUBROUTINE GET_ATTENUATION_PERIODS
          END INTERFACE 
        END MODULE GET_ATTENUATION_PERIODS__genmod
