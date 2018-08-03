        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:35 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_VTK_DATA_POINTS__genmod
          INTERFACE 
            SUBROUTINE WRITE_VTK_DATA_POINTS(NGLOB,XSTORE_DUMMY,        &
     &YSTORE_DUMMY,ZSTORE_DUMMY,POINTS_GLOBALINDICES,                   &
     &NUM_POINTS_GLOBALINDICES,PRNAME_FILE)
              INTEGER(KIND=4) :: NUM_POINTS_GLOBALINDICES
              INTEGER(KIND=4) :: NGLOB
              REAL(KIND=8) :: XSTORE_DUMMY(NGLOB)
              REAL(KIND=8) :: YSTORE_DUMMY(NGLOB)
              REAL(KIND=8) :: ZSTORE_DUMMY(NGLOB)
              INTEGER(KIND=4) :: POINTS_GLOBALINDICES(                  &
     &NUM_POINTS_GLOBALINDICES)
              CHARACTER(LEN=512) :: PRNAME_FILE
            END SUBROUTINE WRITE_VTK_DATA_POINTS
          END INTERFACE 
        END MODULE WRITE_VTK_DATA_POINTS__genmod
