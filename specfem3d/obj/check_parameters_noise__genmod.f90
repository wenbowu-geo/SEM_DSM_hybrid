        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_PARAMETERS_NOISE__genmod
          INTERFACE 
            SUBROUTINE CHECK_PARAMETERS_NOISE(MYRANK,NOISE_TOMOGRAPHY,  &
     &SIMULATION_TYPE,SAVE_FORWARD,LOCAL_PATH,NSPEC_TOP,NSTEP)
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: NOISE_TOMOGRAPHY
              INTEGER(KIND=4) :: SIMULATION_TYPE
              LOGICAL(KIND=4) :: SAVE_FORWARD
              CHARACTER(LEN=512) :: LOCAL_PATH
              INTEGER(KIND=4) :: NSPEC_TOP
              INTEGER(KIND=4) :: NSTEP
            END SUBROUTINE CHECK_PARAMETERS_NOISE
          END INTERFACE 
        END MODULE CHECK_PARAMETERS_NOISE__genmod
