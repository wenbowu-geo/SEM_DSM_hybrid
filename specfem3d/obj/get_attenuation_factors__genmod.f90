        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:32 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ATTENUATION_FACTORS__genmod
          INTERFACE 
            SUBROUTINE GET_ATTENUATION_FACTORS(MYRANK,Q_MU,             &
     &MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,F_C_SOURCE,TAU_SIGMA&
     &,BETA,ONE_MINUS_SUM_BETA,FACTOR_SCALE,Q_KAPPA,BETA_KAPPA,         &
     &ONE_MINUS_SUM_BETA_KAPPA,FACTOR_SCALE_KAPPA,                      &
     &ATTENUATION_F0_REFERENCE)
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: Q_MU
              REAL(KIND=8) :: MIN_ATTENUATION_PERIOD
              REAL(KIND=8) :: MAX_ATTENUATION_PERIOD
              REAL(KIND=8) :: F_C_SOURCE
              REAL(KIND=8) :: TAU_SIGMA(3)
              REAL(KIND=8) :: BETA(3)
              REAL(KIND=8) :: ONE_MINUS_SUM_BETA
              REAL(KIND=8) :: FACTOR_SCALE
              REAL(KIND=8) :: Q_KAPPA
              REAL(KIND=8) :: BETA_KAPPA(3)
              REAL(KIND=8) :: ONE_MINUS_SUM_BETA_KAPPA
              REAL(KIND=8) :: FACTOR_SCALE_KAPPA
              REAL(KIND=8) :: ATTENUATION_F0_REFERENCE
            END SUBROUTINE GET_ATTENUATION_FACTORS
          END INTERFACE 
        END MODULE GET_ATTENUATION_FACTORS__genmod
