        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:58 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SAVE_DEPTH_ID_INNERBOUND__genmod
          INTERFACE 
            SUBROUTINE SAVE_DEPTH_ID_INNERBOUND(MYRANK,NPROC,NELE_BOUND,&
     &ISPEC_IS_ELASTIC,ISPEC_IS_ACOUSTIC,NSPEC,IBOOL,NGLOB,XSTORE_DUMMY,&
     &YSTORE_DUMMY,ZSTORE_DUMMY,MEDIA_TYPE,WZGLL,IMAIN_OUTPUT,PRNAME,   &
     &LOCAL_PATH)
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: NPROC
              INTEGER(KIND=4) :: NELE_BOUND
              LOGICAL(KIND=4) :: ISPEC_IS_ELASTIC(NSPEC)
              LOGICAL(KIND=4) :: ISPEC_IS_ACOUSTIC(NSPEC)
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              REAL(KIND=8) :: XSTORE_DUMMY(NGLOB)
              REAL(KIND=8) :: YSTORE_DUMMY(NGLOB)
              REAL(KIND=8) :: ZSTORE_DUMMY(NGLOB)
              INTEGER(KIND=4) :: MEDIA_TYPE
              REAL(KIND=8) :: WZGLL(5)
              INTEGER(KIND=4) :: IMAIN_OUTPUT
              CHARACTER(LEN=256) :: PRNAME
              CHARACTER(LEN=256) :: LOCAL_PATH
            END SUBROUTINE SAVE_DEPTH_ID_INNERBOUND
          END INTERFACE 
        END MODULE SAVE_DEPTH_ID_INNERBOUND__genmod
