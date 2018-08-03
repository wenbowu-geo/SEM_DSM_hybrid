        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:32 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ATTENUATION_INVERT_BY_SIMPLEX__genmod
          INTERFACE 
            SUBROUTINE ATTENUATION_INVERT_BY_SIMPLEX(T2,T1,N,Q_REAL,    &
     &TAU_S,TAU_EPS)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: T2
              REAL(KIND=8) :: T1
              REAL(KIND=8) :: Q_REAL
              REAL(KIND=8) :: TAU_S(N)
              REAL(KIND=8) :: TAU_EPS(N)
            END SUBROUTINE ATTENUATION_INVERT_BY_SIMPLEX
          END INTERFACE 
        END MODULE ATTENUATION_INVERT_BY_SIMPLEX__genmod
