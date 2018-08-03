        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:28 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MEMORY_EVAL__genmod
          INTERFACE 
            SUBROUTINE MEMORY_EVAL(NSPEC_AB,NGLOB_AB,                   &
     &MAX_NIBOOL_INTERFACES_EXT_MESH,NUM_INTERFACES_EXT_MESH,           &
     &APPROXIMATE_OCEAN_LOAD,MEMORY_SIZE)
              INTEGER(KIND=4), INTENT(IN) :: NSPEC_AB
              INTEGER(KIND=4), INTENT(IN) :: NGLOB_AB
              INTEGER(KIND=4), INTENT(IN) ::                            &
     &MAX_NIBOOL_INTERFACES_EXT_MESH
              INTEGER(KIND=4), INTENT(IN) :: NUM_INTERFACES_EXT_MESH
              LOGICAL(KIND=4), INTENT(IN) :: APPROXIMATE_OCEAN_LOAD
              REAL(KIND=8), INTENT(OUT) :: MEMORY_SIZE
            END SUBROUTINE MEMORY_EVAL
          END INTERFACE 
        END MODULE MEMORY_EVAL__genmod
