        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:17 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_FINAL_PERM__genmod
          INTERFACE 
            SUBROUTINE GET_FINAL_PERM(COLOR,PERM,                       &
     &FIRST_ELEM_NUMBER_IN_THIS_COLOR,NSPEC,NB_COLORS,                  &
     &NB_COLORS_OUTER_ELEMENTS,ISPEC_IS_D,NSPEC_DOMAIN)
              INTEGER(KIND=4), INTENT(IN) :: NB_COLORS
              INTEGER(KIND=4), INTENT(IN) :: NSPEC
              INTEGER(KIND=4), INTENT(IN) :: COLOR(NSPEC)
              INTEGER(KIND=4), INTENT(INOUT) :: PERM(NSPEC)
              INTEGER(KIND=4), INTENT(INOUT) ::                         &
     &FIRST_ELEM_NUMBER_IN_THIS_COLOR(NB_COLORS)
              INTEGER(KIND=4), INTENT(IN) :: NB_COLORS_OUTER_ELEMENTS
              LOGICAL(KIND=4), INTENT(IN) :: ISPEC_IS_D(NSPEC)
              INTEGER(KIND=4), INTENT(IN) :: NSPEC_DOMAIN
            END SUBROUTINE GET_FINAL_PERM
          END INTERFACE 
        END MODULE GET_FINAL_PERM__genmod
