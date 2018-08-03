        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:34 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SAVE_HEADER_FILE__genmod
          INTERFACE 
            SUBROUTINE SAVE_HEADER_FILE(NSPEC_AB,NGLOB_AB,NPROC,        &
     &ATTENUATION,ANISOTROPY,NSTEP,DT,STACEY_INSTEAD_OF_FREE_SURFACE,   &
     &SIMULATION_TYPE,MEMORY_SIZE,NFACES_SURFACE_GLOB_EXT_MESH)
              INTEGER(KIND=4) :: NSPEC_AB
              INTEGER(KIND=4) :: NGLOB_AB
              INTEGER(KIND=4) :: NPROC
              LOGICAL(KIND=4) :: ATTENUATION
              LOGICAL(KIND=4) :: ANISOTROPY
              INTEGER(KIND=4) :: NSTEP
              REAL(KIND=8) :: DT
              LOGICAL(KIND=4) :: STACEY_INSTEAD_OF_FREE_SURFACE
              INTEGER(KIND=4) :: SIMULATION_TYPE
              REAL(KIND=8) :: MEMORY_SIZE
              INTEGER(KIND=4) :: NFACES_SURFACE_GLOB_EXT_MESH
            END SUBROUTINE SAVE_HEADER_FILE
          END INTERFACE 
        END MODULE SAVE_HEADER_FILE__genmod
