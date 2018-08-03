        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:19 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INTERPOL__genmod
          INTERFACE 
            FUNCTION INTERPOL(V,I,X,NL)
              INTEGER(KIND=4) :: NL
              REAL(KIND=8) :: V(NL,4)
              INTEGER(KIND=4) :: I
              REAL(KIND=8) :: X
              REAL(KIND=8) :: INTERPOL
            END FUNCTION INTERPOL
          END INTERFACE 
        END MODULE INTERPOL__genmod
