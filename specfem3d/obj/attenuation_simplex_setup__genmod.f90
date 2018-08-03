        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:32 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ATTENUATION_SIMPLEX_SETUP__genmod
          INTERFACE 
            SUBROUTINE ATTENUATION_SIMPLEX_SETUP(NF_IN,NSLS_IN,F_IN,Q_IN&
     &,TAU_S_IN)
              INTEGER(KIND=4) :: NSLS_IN
              INTEGER(KIND=4) :: NF_IN
              REAL(KIND=8) :: F_IN(NF_IN)
              REAL(KIND=8) :: Q_IN
              REAL(KIND=8) :: TAU_S_IN(NSLS_IN)
            END SUBROUTINE ATTENUATION_SIMPLEX_SETUP
          END INTERFACE 
        END MODULE ATTENUATION_SIMPLEX_SETUP__genmod
