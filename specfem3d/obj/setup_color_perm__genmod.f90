        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:27 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SETUP_COLOR_PERM__genmod
          INTERFACE 
            SUBROUTINE SETUP_COLOR_PERM(MYRANK,NSPEC,NGLOB,IBOOL,       &
     &ANISOTROPY,SAVE_MESH_FILES)
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              LOGICAL(KIND=4) :: ANISOTROPY
              LOGICAL(KIND=4) :: SAVE_MESH_FILES
            END SUBROUTINE SETUP_COLOR_PERM
          END INTERFACE 
        END MODULE SETUP_COLOR_PERM__genmod
