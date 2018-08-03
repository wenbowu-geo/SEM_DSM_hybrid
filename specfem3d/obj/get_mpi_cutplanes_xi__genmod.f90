        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:48 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_MPI_CUTPLANES_XI__genmod
          INTERFACE 
            SUBROUTINE GET_MPI_CUTPLANES_XI(MYRANK,PRNAME,NSPEC,        &
     &IMPICUT_XI,IBOOL,XSTORE,YSTORE,ZSTORE,MASK_IBOOL,NPOINTOT,        &
     &NSPEC2D_A_ETA,NSPEC2D_B_ETA)
              INTEGER(KIND=4) :: NPOINTOT
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              CHARACTER(LEN=512) :: PRNAME
              LOGICAL(KIND=4) :: IMPICUT_XI(2,NSPEC)
              INTEGER(KIND=4) :: IBOOL(2,2,2,NSPEC)
              REAL(KIND=8) :: XSTORE(2,2,2,NSPEC)
              REAL(KIND=8) :: YSTORE(2,2,2,NSPEC)
              REAL(KIND=8) :: ZSTORE(2,2,2,NSPEC)
              LOGICAL(KIND=4) :: MASK_IBOOL(NPOINTOT)
              INTEGER(KIND=4) :: NSPEC2D_A_ETA
              INTEGER(KIND=4) :: NSPEC2D_B_ETA
            END SUBROUTINE GET_MPI_CUTPLANES_XI
          END INTERFACE 
        END MODULE GET_MPI_CUTPLANES_XI__genmod
