        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:35 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_VTK_DATA_ELEM_VECTORS__genmod
          INTERFACE 
            SUBROUTINE WRITE_VTK_DATA_ELEM_VECTORS(NSPEC,NGLOB,         &
     &XSTORE_DUMMY,YSTORE_DUMMY,ZSTORE_DUMMY,IBOOL,ELEM_VECTOR,         &
     &PRNAME_FILE)
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: NSPEC
              REAL(KIND=8) :: XSTORE_DUMMY(NGLOB)
              REAL(KIND=8) :: YSTORE_DUMMY(NGLOB)
              REAL(KIND=8) :: ZSTORE_DUMMY(NGLOB)
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              REAL(KIND=8) :: ELEM_VECTOR(3,NSPEC)
              CHARACTER(LEN=512) :: PRNAME_FILE
            END SUBROUTINE WRITE_VTK_DATA_ELEM_VECTORS
          END INTERFACE 
        END MODULE WRITE_VTK_DATA_ELEM_VECTORS__genmod
