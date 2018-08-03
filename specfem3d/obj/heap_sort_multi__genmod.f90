        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:38 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HEAP_SORT_MULTI__genmod
          INTERFACE 
            SUBROUTINE HEAP_SORT_MULTI(N,DX,DY,DZ,IA,IB)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(INOUT) :: DX(N)
              REAL(KIND=8), INTENT(INOUT) :: DY(N)
              REAL(KIND=8), INTENT(INOUT) :: DZ(N)
              INTEGER(KIND=4), INTENT(INOUT) :: IA(N)
              INTEGER(KIND=4), INTENT(INOUT) :: IB(N)
            END SUBROUTINE HEAP_SORT_MULTI
          END INTERFACE 
        END MODULE HEAP_SORT_MULTI__genmod
