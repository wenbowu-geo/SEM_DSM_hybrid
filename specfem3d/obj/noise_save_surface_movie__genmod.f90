        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NOISE_SAVE_SURFACE_MOVIE__genmod
          INTERFACE 
            SUBROUTINE NOISE_SAVE_SURFACE_MOVIE(DISPL,IBOOL,            &
     &NOISE_SURFACE_MOVIE,IT,NSPEC_AB_VAL,NGLOB_AB_VAL,                 &
     &NUM_FREE_SURFACE_FACES,FREE_SURFACE_ISPEC,FREE_SURFACE_IJK,       &
     &MESH_POINTER,GPU_MODE)
              INTEGER(KIND=4) :: NUM_FREE_SURFACE_FACES
              INTEGER(KIND=4) :: NGLOB_AB_VAL
              INTEGER(KIND=4) :: NSPEC_AB_VAL
              REAL(KIND=8) :: DISPL(3,NGLOB_AB_VAL)
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC_AB_VAL)
              REAL(KIND=8) :: NOISE_SURFACE_MOVIE(3,25,                 &
     &NUM_FREE_SURFACE_FACES)
              INTEGER(KIND=4) :: IT
              INTEGER(KIND=4) :: FREE_SURFACE_ISPEC(                    &
     &NUM_FREE_SURFACE_FACES)
              INTEGER(KIND=4) :: FREE_SURFACE_IJK(3,25,                 &
     &NUM_FREE_SURFACE_FACES)
              INTEGER(KIND=8) :: MESH_POINTER
              LOGICAL(KIND=4) :: GPU_MODE
            END SUBROUTINE NOISE_SAVE_SURFACE_MOVIE
          END INTERFACE 
        END MODULE NOISE_SAVE_SURFACE_MOVIE__genmod
