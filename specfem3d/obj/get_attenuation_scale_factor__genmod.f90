        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:32 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ATTENUATION_SCALE_FACTOR__genmod
          INTERFACE 
            SUBROUTINE GET_ATTENUATION_SCALE_FACTOR(MYRANK,F_C_SOURCE,  &
     &TAU_EPS,TAU_SIGMA,Q_VAL,SCALE_FACTOR,ATTENUATION_F0_REFERENCE)
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: F_C_SOURCE
              REAL(KIND=8) :: TAU_EPS(3)
              REAL(KIND=8) :: TAU_SIGMA(3)
              REAL(KIND=8) :: Q_VAL
              REAL(KIND=8) :: SCALE_FACTOR
              REAL(KIND=8), INTENT(IN) :: ATTENUATION_F0_REFERENCE
            END SUBROUTINE GET_ATTENUATION_SCALE_FACTOR
          END INTERFACE 
        END MODULE GET_ATTENUATION_SCALE_FACTOR__genmod
