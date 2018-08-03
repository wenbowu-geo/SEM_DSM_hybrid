        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:28 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_MESH_DISTANCES__genmod
          INTERFACE 
            SUBROUTINE CHECK_MESH_DISTANCES(MYRANK,NSPEC_AB,NGLOB_AB,   &
     &IBOOL,XSTORE,YSTORE,ZSTORE,X_MIN_GLOB,X_MAX_GLOB,Y_MIN_GLOB,      &
     &Y_MAX_GLOB,Z_MIN_GLOB,Z_MAX_GLOB,ELEMSIZE_MIN_GLOB,               &
     &ELEMSIZE_MAX_GLOB,DISTANCE_MIN_GLOB,DISTANCE_MAX_GLOB)
              INTEGER(KIND=4), INTENT(IN) :: NGLOB_AB
              INTEGER(KIND=4), INTENT(IN) :: NSPEC_AB
              INTEGER(KIND=4), INTENT(IN) :: MYRANK
              INTEGER(KIND=4), INTENT(IN) :: IBOOL(5,5,5,NSPEC_AB)
              REAL(KIND=8), INTENT(IN) :: XSTORE(NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: YSTORE(NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: ZSTORE(NGLOB_AB)
              REAL(KIND=8), INTENT(OUT) :: X_MIN_GLOB
              REAL(KIND=8), INTENT(OUT) :: X_MAX_GLOB
              REAL(KIND=8), INTENT(OUT) :: Y_MIN_GLOB
              REAL(KIND=8), INTENT(OUT) :: Y_MAX_GLOB
              REAL(KIND=8), INTENT(OUT) :: Z_MIN_GLOB
              REAL(KIND=8), INTENT(OUT) :: Z_MAX_GLOB
              REAL(KIND=8), INTENT(OUT) :: ELEMSIZE_MIN_GLOB
              REAL(KIND=8), INTENT(OUT) :: ELEMSIZE_MAX_GLOB
              REAL(KIND=8), INTENT(OUT) :: DISTANCE_MIN_GLOB
              REAL(KIND=8), INTENT(OUT) :: DISTANCE_MAX_GLOB
            END SUBROUTINE CHECK_MESH_DISTANCES
          END INTERFACE 
        END MODULE CHECK_MESH_DISTANCES__genmod
