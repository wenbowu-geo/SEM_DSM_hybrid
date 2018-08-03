        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SAVE_KERNELS_STRENGTH_NOISE__genmod
          INTERFACE 
            SUBROUTINE SAVE_KERNELS_STRENGTH_NOISE(MYRANK,LOCAL_PATH,   &
     &SIGMA_KL,NSPEC_AB_VAL)
              INTEGER(KIND=4) :: NSPEC_AB_VAL
              INTEGER(KIND=4) :: MYRANK
              CHARACTER(LEN=512) :: LOCAL_PATH
              REAL(KIND=8) :: SIGMA_KL(5,5,5,NSPEC_AB_VAL)
            END SUBROUTINE SAVE_KERNELS_STRENGTH_NOISE
          END INTERFACE 
        END MODULE SAVE_KERNELS_STRENGTH_NOISE__genmod
