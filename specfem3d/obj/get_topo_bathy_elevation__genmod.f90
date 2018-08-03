        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:52 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_TOPO_BATHY_ELEVATION__genmod
          INTERFACE 
            SUBROUTINE GET_TOPO_BATHY_ELEVATION(X_TARGET,Y_TARGET,      &
     &TARGET_ELEVATION,ITOPO_BATHY,NX_TOPO,NY_TOPO,UTM_PROJECTION_ZONE, &
     &SUPPRESS_UTM_PROJECTION)
              INTEGER(KIND=4) :: NY_TOPO
              INTEGER(KIND=4) :: NX_TOPO
              REAL(KIND=8), INTENT(IN) :: X_TARGET
              REAL(KIND=8), INTENT(IN) :: Y_TARGET
              REAL(KIND=8), INTENT(OUT) :: TARGET_ELEVATION
              INTEGER(KIND=4) :: ITOPO_BATHY(NX_TOPO,NY_TOPO)
              INTEGER(KIND=4) :: UTM_PROJECTION_ZONE
              LOGICAL(KIND=4) :: SUPPRESS_UTM_PROJECTION
            END SUBROUTINE GET_TOPO_BATHY_ELEVATION
          END INTERFACE 
        END MODULE GET_TOPO_BATHY_ELEVATION__genmod
