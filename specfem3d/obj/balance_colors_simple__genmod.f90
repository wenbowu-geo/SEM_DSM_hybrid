        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:17 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BALANCE_COLORS_SIMPLE__genmod
          INTERFACE 
            SUBROUTINE BALANCE_COLORS_SIMPLE(IBOOL,IS_ON_A_SLICE_EDGE,  &
     &ISPEC_IS_D,MYRANK,NSPEC,NGLOB,COLOR,NB_COLORS_OUTER_ELEMENTS,     &
     &NB_COLORS_INNER_ELEMENTS,NSPEC_OUTER,NSPEC_INNER,MASK_IBOOL)
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              LOGICAL(KIND=4) :: IS_ON_A_SLICE_EDGE(NSPEC)
              LOGICAL(KIND=4) :: ISPEC_IS_D(NSPEC)
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: COLOR(NSPEC)
              INTEGER(KIND=4) :: NB_COLORS_OUTER_ELEMENTS
              INTEGER(KIND=4) :: NB_COLORS_INNER_ELEMENTS
              INTEGER(KIND=4) :: NSPEC_OUTER
              INTEGER(KIND=4) :: NSPEC_INNER
              LOGICAL(KIND=4) :: MASK_IBOOL(NGLOB)
            END SUBROUTINE BALANCE_COLORS_SIMPLE
          END INTERFACE 
        END MODULE BALANCE_COLORS_SIMPLE__genmod
