        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:45 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CVD_WRITE_CORNER_ELEMENTS__genmod
          INTERFACE 
            SUBROUTINE CVD_WRITE_CORNER_ELEMENTS(NSPEC_AB,NGLOB_AB,IBOOL&
     &,NP,NELEMENT,IT,NEE,NUMPOIN)
              INTEGER(KIND=4), INTENT(IN) :: NSPEC_AB
              INTEGER(KIND=4), INTENT(IN) :: NGLOB_AB
              INTEGER(KIND=4), INTENT(IN) :: IBOOL(5,5,5,NSPEC_AB)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: NELEMENT
              INTEGER(KIND=4) :: IT
              INTEGER(KIND=4) :: NEE
              INTEGER(KIND=4) :: NUMPOIN
            END SUBROUTINE CVD_WRITE_CORNER_ELEMENTS
          END INTERFACE 
        END MODULE CVD_WRITE_CORNER_ELEMENTS__genmod
