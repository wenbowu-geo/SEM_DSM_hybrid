        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:45:17 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INVERT_HESS__genmod
          INTERFACE 
            SUBROUTINE INVERT_HESS(HESS_MATRIX)
              USE TOMOGRAPHY_PAR
              REAL(KIND=8) :: HESS_MATRIX(5,5,5,NSPEC)
            END SUBROUTINE INVERT_HESS
          END INTERFACE 
        END MODULE INVERT_HESS__genmod
