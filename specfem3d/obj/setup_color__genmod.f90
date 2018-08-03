        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:27 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SETUP_COLOR__genmod
          INTERFACE 
            SUBROUTINE SETUP_COLOR(MYRANK,NSPEC,NGLOB,IBOOL,PERM,       &
     &ISPEC_IS_D,IDOMAIN,NUM_PHASE_ISPEC_D,PHASE_ISPEC_INNER_D,         &
     &SAVE_MESH_FILES)
              INTEGER(KIND=4) :: NUM_PHASE_ISPEC_D
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              INTEGER(KIND=4) :: PERM(NSPEC)
              LOGICAL(KIND=4) :: ISPEC_IS_D(NSPEC)
              INTEGER(KIND=4) :: IDOMAIN
              INTEGER(KIND=4) :: PHASE_ISPEC_INNER_D(NUM_PHASE_ISPEC_D,2&
     &)
              LOGICAL(KIND=4) :: SAVE_MESH_FILES
            END SUBROUTINE SETUP_COLOR
          END INTERFACE 
        END MODULE SETUP_COLOR__genmod
