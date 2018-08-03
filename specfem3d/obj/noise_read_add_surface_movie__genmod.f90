        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NOISE_READ_ADD_SURFACE_MOVIE__genmod
          INTERFACE 
            SUBROUTINE NOISE_READ_ADD_SURFACE_MOVIE(NMOVIE_POINTS,ACCEL,&
     &NORMAL_X_NOISE,NORMAL_Y_NOISE,NORMAL_Z_NOISE,MASK_NOISE,IBOOL,    &
     &NOISE_SURFACE_MOVIE,IT,NSPEC_AB_VAL,NGLOB_AB_VAL,                 &
     &NUM_FREE_SURFACE_FACES,FREE_SURFACE_ISPEC,FREE_SURFACE_IJK,       &
     &FREE_SURFACE_JACOBIAN2DW)
              INTEGER(KIND=4) :: NUM_FREE_SURFACE_FACES
              INTEGER(KIND=4) :: NGLOB_AB_VAL
              INTEGER(KIND=4) :: NSPEC_AB_VAL
              INTEGER(KIND=4) :: NMOVIE_POINTS
              REAL(KIND=8) :: ACCEL(3,NGLOB_AB_VAL)
              REAL(KIND=8) :: NORMAL_X_NOISE(NMOVIE_POINTS)
              REAL(KIND=8) :: NORMAL_Y_NOISE(NMOVIE_POINTS)
              REAL(KIND=8) :: NORMAL_Z_NOISE(NMOVIE_POINTS)
              REAL(KIND=8) :: MASK_NOISE(NMOVIE_POINTS)
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC_AB_VAL)
              REAL(KIND=8) :: NOISE_SURFACE_MOVIE(3,25,                 &
     &NUM_FREE_SURFACE_FACES)
              INTEGER(KIND=4) :: IT
              INTEGER(KIND=4) :: FREE_SURFACE_ISPEC(                    &
     &NUM_FREE_SURFACE_FACES)
              INTEGER(KIND=4) :: FREE_SURFACE_IJK(3,25,                 &
     &NUM_FREE_SURFACE_FACES)
              REAL(KIND=8) :: FREE_SURFACE_JACOBIAN2DW(25,              &
     &NUM_FREE_SURFACE_FACES)
            END SUBROUTINE NOISE_READ_ADD_SURFACE_MOVIE
          END INTERFACE 
        END MODULE NOISE_READ_ADD_SURFACE_MOVIE__genmod
