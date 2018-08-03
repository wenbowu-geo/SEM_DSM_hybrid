        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:32 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ATTENUATION_TAU_SIGMA__genmod
          INTERFACE 
            SUBROUTINE GET_ATTENUATION_TAU_SIGMA(TAU_S,NSLS,MIN_PERIOD, &
     &MAX_PERIOD)
              INTEGER(KIND=4) :: NSLS
              REAL(KIND=8), INTENT(OUT) :: TAU_S(NSLS)
              REAL(KIND=8), INTENT(IN) :: MIN_PERIOD
              REAL(KIND=8), INTENT(IN) :: MAX_PERIOD
            END SUBROUTINE GET_ATTENUATION_TAU_SIGMA
          END INTERFACE 
        END MODULE GET_ATTENUATION_TAU_SIGMA__genmod
