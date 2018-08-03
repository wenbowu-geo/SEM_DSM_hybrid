        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INTGRL__genmod
          INTERFACE 
            SUBROUTINE INTGRL(SUMVAL,R,NIR,NER,F,S1,S2,S3)
              REAL(KIND=8) :: SUMVAL
              REAL(KIND=8) :: R(640)
              INTEGER(KIND=4) :: NIR
              INTEGER(KIND=4) :: NER
              REAL(KIND=8) :: F(640)
              REAL(KIND=8) :: S1(640)
              REAL(KIND=8) :: S2(640)
              REAL(KIND=8) :: S3(640)
            END SUBROUTINE INTGRL
          END INTERFACE 
        END MODULE INTGRL__genmod
