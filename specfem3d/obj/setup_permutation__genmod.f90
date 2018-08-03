        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:27 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SETUP_PERMUTATION__genmod
          INTERFACE 
            SUBROUTINE SETUP_PERMUTATION(MYRANK,NSPEC,NGLOB,IBOOL,      &
     &ANISOTROPY,PERM,SAVE_MESH_FILES)
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              LOGICAL(KIND=4) :: ANISOTROPY
              INTEGER(KIND=4), INTENT(INOUT) :: PERM(NSPEC)
              LOGICAL(KIND=4) :: SAVE_MESH_FILES
            END SUBROUTINE SETUP_PERMUTATION
          END INTERFACE 
        END MODULE SETUP_PERMUTATION__genmod
