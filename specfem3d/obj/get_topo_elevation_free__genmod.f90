        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:52 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_TOPO_ELEVATION_FREE__genmod
          INTERFACE 
            SUBROUTINE GET_TOPO_ELEVATION_FREE(X_TARGET,Y_TARGET,       &
     &TARGET_ELEVATION,TARGET_DISTMIN,NSPEC_AB,NGLOB_AB,IBOOL,XSTORE,   &
     &YSTORE,ZSTORE,NUM_FREE_SURFACE_FACES,FREE_SURFACE_ISPEC,          &
     &FREE_SURFACE_IJK)
              INTEGER(KIND=4) :: NUM_FREE_SURFACE_FACES
              INTEGER(KIND=4) :: NGLOB_AB
              INTEGER(KIND=4) :: NSPEC_AB
              REAL(KIND=8), INTENT(IN) :: X_TARGET
              REAL(KIND=8), INTENT(IN) :: Y_TARGET
              REAL(KIND=8), INTENT(OUT) :: TARGET_ELEVATION
              REAL(KIND=8), INTENT(OUT) :: TARGET_DISTMIN
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: XSTORE(NGLOB_AB)
              REAL(KIND=8) :: YSTORE(NGLOB_AB)
              REAL(KIND=8) :: ZSTORE(NGLOB_AB)
              INTEGER(KIND=4) :: FREE_SURFACE_ISPEC(                    &
     &NUM_FREE_SURFACE_FACES)
              INTEGER(KIND=4) :: FREE_SURFACE_IJK(3,25,                 &
     &NUM_FREE_SURFACE_FACES)
            END SUBROUTINE GET_TOPO_ELEVATION_FREE
          END INTERFACE 
        END MODULE GET_TOPO_ELEVATION_FREE__genmod
