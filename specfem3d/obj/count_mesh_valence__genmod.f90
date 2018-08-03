        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:17 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COUNT_MESH_VALENCE__genmod
          INTERFACE 
            SUBROUTINE COUNT_MESH_VALENCE(IBOOL,IS_ON_A_SLICE_EDGE,     &
     &ISPEC_IS_D,MYRANK,NSPEC,NGLOB,MAXVAL_COUNT_IBOOL_OUTER,           &
     &MAXVAL_COUNT_IBOOL_INNER)
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              LOGICAL(KIND=4) :: IS_ON_A_SLICE_EDGE(NSPEC)
              LOGICAL(KIND=4) :: ISPEC_IS_D(NSPEC)
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: MAXVAL_COUNT_IBOOL_OUTER
              INTEGER(KIND=4) :: MAXVAL_COUNT_IBOOL_INNER
            END SUBROUTINE COUNT_MESH_VALENCE
          END INTERFACE 
        END MODULE COUNT_MESH_VALENCE__genmod
