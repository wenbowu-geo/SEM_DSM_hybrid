        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:59 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_STACEY_ACOUSTIC_GPU__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_STACEY_ACOUSTIC_GPU(IPHASE,              &
     &NUM_ABS_BOUNDARY_FACES,SIMULATION_TYPE,SAVE_FORWARD,NSTEP,IT,     &
     &B_RECLEN_POTENTIAL,B_ABSORB_POTENTIAL,B_NUM_ABS_BOUNDARY_FACES,   &
     &MESH_POINTER)
              INTEGER(KIND=4) :: B_NUM_ABS_BOUNDARY_FACES
              INTEGER(KIND=4) :: IPHASE
              INTEGER(KIND=4) :: NUM_ABS_BOUNDARY_FACES
              INTEGER(KIND=4) :: SIMULATION_TYPE
              LOGICAL(KIND=4) :: SAVE_FORWARD
              INTEGER(KIND=4) :: NSTEP
              INTEGER(KIND=4) :: IT
              INTEGER(KIND=4) :: B_RECLEN_POTENTIAL
              REAL(KIND=8) :: B_ABSORB_POTENTIAL(25,                    &
     &B_NUM_ABS_BOUNDARY_FACES)
              INTEGER(KIND=8) :: MESH_POINTER
            END SUBROUTINE COMPUTE_STACEY_ACOUSTIC_GPU
          END INTERFACE 
        END MODULE COMPUTE_STACEY_ACOUSTIC_GPU__genmod
