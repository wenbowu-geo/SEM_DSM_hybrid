        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:32 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ATTENUATION_MAXWELL__genmod
          INTERFACE 
            SUBROUTINE ATTENUATION_MAXWELL(NF,NSLS,F,TAU_S,TAU_EPS,B,A)
              INTEGER(KIND=4) :: NSLS
              INTEGER(KIND=4) :: NF
              REAL(KIND=8) :: F(NF)
              REAL(KIND=8) :: TAU_S(NSLS)
              REAL(KIND=8) :: TAU_EPS(NSLS)
              REAL(KIND=8) :: B(NF)
              REAL(KIND=8) :: A(NF)
            END SUBROUTINE ATTENUATION_MAXWELL
          END INTERFACE 
        END MODULE ATTENUATION_MAXWELL__genmod
