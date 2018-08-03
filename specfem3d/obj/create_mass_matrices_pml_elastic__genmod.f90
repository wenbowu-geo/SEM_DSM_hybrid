        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:54 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CREATE_MASS_MATRICES_PML_ELASTIC__genmod
          INTERFACE 
            SUBROUTINE CREATE_MASS_MATRICES_PML_ELASTIC(NSPEC,IBOOL)
              INTEGER(KIND=4), INTENT(IN) :: NSPEC
              INTEGER(KIND=4), INTENT(IN) :: IBOOL(5,5,5,NSPEC)
            END SUBROUTINE CREATE_MASS_MATRICES_PML_ELASTIC
          END INTERFACE 
        END MODULE CREATE_MASS_MATRICES_PML_ELASTIC__genmod
