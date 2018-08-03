        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ADD_SOURCE_MASTER_REC_NOISE__genmod
          INTERFACE 
            SUBROUTINE ADD_SOURCE_MASTER_REC_NOISE(MYRANK,NREC,NSTEP,   &
     &ACCEL,NOISE_SOURCEARRAY,IBOOL,ISLICE_SELECTED_REC,                &
     &ISPEC_SELECTED_REC,IT,IREC_MASTER_NOISE,NSPEC_AB_VAL,NGLOB_AB_VAL)
              INTEGER(KIND=4) :: NGLOB_AB_VAL
              INTEGER(KIND=4) :: NSPEC_AB_VAL
              INTEGER(KIND=4) :: NSTEP
              INTEGER(KIND=4) :: NREC
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: ACCEL(3,NGLOB_AB_VAL)
              REAL(KIND=8) :: NOISE_SOURCEARRAY(3,5,5,5,NSTEP)
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC_AB_VAL)
              INTEGER(KIND=4) :: ISLICE_SELECTED_REC(NREC)
              INTEGER(KIND=4) :: ISPEC_SELECTED_REC(NREC)
              INTEGER(KIND=4) :: IT
              INTEGER(KIND=4) :: IREC_MASTER_NOISE
            END SUBROUTINE ADD_SOURCE_MASTER_REC_NOISE
          END INTERFACE 
        END MODULE ADD_SOURCE_MASTER_REC_NOISE__genmod
