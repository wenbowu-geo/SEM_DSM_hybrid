        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:00 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_STACEY_VISCOELASTIC_GPU__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_STACEY_VISCOELASTIC_GPU(IPHASE,          &
     &NUM_ABS_BOUNDARY_FACES,SIMULATION_TYPE,SAVE_FORWARD,NSTEP,IT,     &
     &B_NUM_ABS_BOUNDARY_FACES,B_RECLEN_FIELD,B_ABSORB_FIELD,           &
     &MESH_POINTER)
              INTEGER(KIND=4) :: B_NUM_ABS_BOUNDARY_FACES
              INTEGER(KIND=4) :: IPHASE
              INTEGER(KIND=4) :: NUM_ABS_BOUNDARY_FACES
              INTEGER(KIND=4) :: SIMULATION_TYPE
              LOGICAL(KIND=4) :: SAVE_FORWARD
              INTEGER(KIND=4) :: NSTEP
              INTEGER(KIND=4) :: IT
              INTEGER(KIND=4) :: B_RECLEN_FIELD
              REAL(KIND=8) :: B_ABSORB_FIELD(3,25,                      &
     &B_NUM_ABS_BOUNDARY_FACES)
              INTEGER(KIND=8) :: MESH_POINTER
            END SUBROUTINE COMPUTE_STACEY_VISCOELASTIC_GPU
          END INTERFACE 
        END MODULE COMPUTE_STACEY_VISCOELASTIC_GPU__genmod
