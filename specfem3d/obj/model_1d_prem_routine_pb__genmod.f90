        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MODEL_1D_PREM_ROUTINE_PB__genmod
          INTERFACE 
            SUBROUTINE MODEL_1D_PREM_ROUTINE_PB(XLOC,YLOC,ZLOC,RO_PREM, &
     &VP_PREM,VS_PREM,IDOM)
              REAL(KIND=8), INTENT(IN) :: XLOC
              REAL(KIND=8), INTENT(IN) :: YLOC
              REAL(KIND=8), INTENT(IN) :: ZLOC
              REAL(KIND=8), INTENT(INOUT) :: RO_PREM
              REAL(KIND=8), INTENT(INOUT) :: VP_PREM
              REAL(KIND=8), INTENT(INOUT) :: VS_PREM
              INTEGER(KIND=4), INTENT(IN) :: IDOM
            END SUBROUTINE MODEL_1D_PREM_ROUTINE_PB
          END INTERFACE 
        END MODULE MODEL_1D_PREM_ROUTINE_PB__genmod
