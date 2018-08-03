        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:52 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UTM_GEO__genmod
          INTERFACE 
            SUBROUTINE UTM_GEO(RLON4,RLAT4,RX4,RY4,UTM_PROJECTION_ZONE, &
     &IWAY,SUPPRESS_UTM_PROJECTION)
              REAL(KIND=8) :: RLON4
              REAL(KIND=8) :: RLAT4
              REAL(KIND=8) :: RX4
              REAL(KIND=8) :: RY4
              INTEGER(KIND=4), INTENT(IN) :: UTM_PROJECTION_ZONE
              INTEGER(KIND=4), INTENT(IN) :: IWAY
              LOGICAL(KIND=4), INTENT(IN) :: SUPPRESS_UTM_PROJECTION
            END SUBROUTINE UTM_GEO
          END INTERFACE 
        END MODULE UTM_GEO__genmod
