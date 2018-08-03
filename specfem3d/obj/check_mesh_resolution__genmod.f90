        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:28 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_MESH_RESOLUTION__genmod
          INTERFACE 
            SUBROUTINE CHECK_MESH_RESOLUTION(MYRANK,NSPEC_AB,NGLOB_AB,  &
     &IBOOL,XSTORE,YSTORE,ZSTORE,KAPPASTORE,MUSTORE,RHO_VP,RHO_VS,DT,   &
     &MODEL_SPEED_MAX,MIN_RESOLVED_PERIOD,LOCAL_PATH,SAVE_MESH_FILES)
              INTEGER(KIND=4) :: NGLOB_AB
              INTEGER(KIND=4) :: NSPEC_AB
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: XSTORE(NGLOB_AB)
              REAL(KIND=8) :: YSTORE(NGLOB_AB)
              REAL(KIND=8) :: ZSTORE(NGLOB_AB)
              REAL(KIND=8) :: KAPPASTORE(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: MUSTORE(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHO_VP(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHO_VS(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: MODEL_SPEED_MAX
              REAL(KIND=8) :: MIN_RESOLVED_PERIOD
              CHARACTER(LEN=512) :: LOCAL_PATH
              LOGICAL(KIND=4) :: SAVE_MESH_FILES
            END SUBROUTINE CHECK_MESH_RESOLUTION
          END INTERFACE 
        END MODULE CHECK_MESH_RESOLUTION__genmod
