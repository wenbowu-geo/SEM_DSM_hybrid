        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NOISE_READ_ADD_SURFACE_MOVIE_GPU__genmod
          INTERFACE 
            SUBROUTINE NOISE_READ_ADD_SURFACE_MOVIE_GPU(                &
     &NOISE_SURFACE_MOVIE,IT,NUM_FREE_SURFACE_FACES,MESH_POINTER,       &
     &NOISE_TOMOGRAPHY)
              INTEGER(KIND=4) :: NUM_FREE_SURFACE_FACES
              REAL(KIND=8) :: NOISE_SURFACE_MOVIE(3,25,                 &
     &NUM_FREE_SURFACE_FACES)
              INTEGER(KIND=4) :: IT
              INTEGER(KIND=8) :: MESH_POINTER
              INTEGER(KIND=4) :: NOISE_TOMOGRAPHY
            END SUBROUTINE NOISE_READ_ADD_SURFACE_MOVIE_GPU
          END INTERFACE 
        END MODULE NOISE_READ_ADD_SURFACE_MOVIE_GPU__genmod
