        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:28 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_MESH_RESOLUTION_PORO__genmod
          INTERFACE 
            SUBROUTINE CHECK_MESH_RESOLUTION_PORO(MYRANK,NSPEC_AB,      &
     &NGLOB_AB,IBOOL,XSTORE,YSTORE,ZSTORE,DT,MODEL_SPEED_MAX,           &
     &MIN_RESOLVED_PERIOD,PHISTORE,TORTSTORE,RHOARRAYSTORE,RHO_VPI,     &
     &RHO_VPII,RHO_VSI,LOCAL_PATH,SAVE_MESH_FILES)
              INTEGER(KIND=4) :: NGLOB_AB
              INTEGER(KIND=4) :: NSPEC_AB
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: XSTORE(NGLOB_AB)
              REAL(KIND=8) :: YSTORE(NGLOB_AB)
              REAL(KIND=8) :: ZSTORE(NGLOB_AB)
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: MODEL_SPEED_MAX
              REAL(KIND=8) :: MIN_RESOLVED_PERIOD
              REAL(KIND=8) :: PHISTORE(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: TORTSTORE(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHOARRAYSTORE(2,5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHO_VPI(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHO_VPII(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHO_VSI(5,5,5,NSPEC_AB)
              CHARACTER(LEN=512) :: LOCAL_PATH
              LOGICAL(KIND=4) :: SAVE_MESH_FILES
            END SUBROUTINE CHECK_MESH_RESOLUTION_PORO
          END INTERFACE 
        END MODULE CHECK_MESH_RESOLUTION_PORO__genmod
