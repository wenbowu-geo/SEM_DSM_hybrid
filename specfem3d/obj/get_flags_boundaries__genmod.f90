        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:48 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_FLAGS_BOUNDARIES__genmod
          INTERFACE 
            SUBROUTINE GET_FLAGS_BOUNDARIES(NSPEC,IPROC_XI,IPROC_ETA,   &
     &ISPEC,IDOUBLING,XSTORE,YSTORE,ZSTORE,IBOUN,IMPICUT_XI,IMPICUT_ETA,&
     &NPROC_XI,NPROC_ETA,UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,       &
     &Z_DEPTH_BLOCK,NEX_XI,NEX_ETA)
              INTEGER(KIND=4), INTENT(IN) :: NSPEC
              INTEGER(KIND=4), INTENT(IN) :: IPROC_XI
              INTEGER(KIND=4), INTENT(IN) :: IPROC_ETA
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
              INTEGER(KIND=4), INTENT(IN) :: IDOUBLING
              REAL(KIND=8), INTENT(IN) :: XSTORE(2,2,2)
              REAL(KIND=8), INTENT(IN) :: YSTORE(2,2,2)
              REAL(KIND=8), INTENT(IN) :: ZSTORE(2,2,2)
              LOGICAL(KIND=4), INTENT(INOUT) :: IBOUN(6,NSPEC)
              LOGICAL(KIND=4), INTENT(INOUT) :: IMPICUT_XI(2,NSPEC)
              LOGICAL(KIND=4), INTENT(INOUT) :: IMPICUT_ETA(2,NSPEC)
              INTEGER(KIND=4), INTENT(IN) :: NPROC_XI
              INTEGER(KIND=4), INTENT(IN) :: NPROC_ETA
              REAL(KIND=8), INTENT(IN) :: UTM_X_MIN
              REAL(KIND=8), INTENT(IN) :: UTM_X_MAX
              REAL(KIND=8), INTENT(IN) :: UTM_Y_MIN
              REAL(KIND=8), INTENT(IN) :: UTM_Y_MAX
              REAL(KIND=8), INTENT(IN) :: Z_DEPTH_BLOCK
              INTEGER(KIND=4), INTENT(IN) :: NEX_XI
              INTEGER(KIND=4), INTENT(IN) :: NEX_ETA
            END SUBROUTINE GET_FLAGS_BOUNDARIES
          END INTERFACE 
        END MODULE GET_FLAGS_BOUNDARIES__genmod
