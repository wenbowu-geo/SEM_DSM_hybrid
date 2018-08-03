        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:37 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_ELASTIC_LDDRK__genmod
          INTERFACE 
            SUBROUTINE UPDATE_ELASTIC_LDDRK(NGLOB,NGLOB_LDDRK,DISPL,    &
     &VELOC,ACCEL,DISPL_LDDRK,VELOC_LDDRK,DELTAT,ALPHA,BETA)
              INTEGER(KIND=4), INTENT(IN) :: NGLOB_LDDRK
              INTEGER(KIND=4), INTENT(IN) :: NGLOB
              REAL(KIND=8), INTENT(INOUT) :: DISPL(3,NGLOB)
              REAL(KIND=8), INTENT(INOUT) :: VELOC(3,NGLOB)
              REAL(KIND=8), INTENT(INOUT) :: ACCEL(3,NGLOB)
              REAL(KIND=8), INTENT(INOUT) :: DISPL_LDDRK(3,NGLOB_LDDRK)
              REAL(KIND=8), INTENT(INOUT) :: VELOC_LDDRK(3,NGLOB_LDDRK)
              REAL(KIND=8), INTENT(IN) :: DELTAT
              REAL(KIND=8), INTENT(IN) :: ALPHA
              REAL(KIND=8), INTENT(IN) :: BETA
            END SUBROUTINE UPDATE_ELASTIC_LDDRK
          END INTERFACE 
        END MODULE UPDATE_ELASTIC_LDDRK__genmod
