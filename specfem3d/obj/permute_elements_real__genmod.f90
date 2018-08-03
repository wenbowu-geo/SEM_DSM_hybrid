        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:17 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PERMUTE_ELEMENTS_REAL__genmod
          INTERFACE 
            SUBROUTINE PERMUTE_ELEMENTS_REAL(ARRAY_TO_PERMUTE,TEMP_ARRAY&
     &,PERM,NSPEC)
              INTEGER(KIND=4), INTENT(IN) :: NSPEC
              REAL(KIND=8), INTENT(INOUT) :: ARRAY_TO_PERMUTE(5,5,5,    &
     &NSPEC)
              REAL(KIND=8), INTENT(INOUT) :: TEMP_ARRAY(5,5,5,NSPEC)
              INTEGER(KIND=4), INTENT(IN) :: PERM(NSPEC)
            END SUBROUTINE PERMUTE_ELEMENTS_REAL
          END INTERFACE 
        END MODULE PERMUTE_ELEMENTS_REAL__genmod
