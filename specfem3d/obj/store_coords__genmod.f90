        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:50 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STORE_COORDS__genmod
          INTERFACE 
            SUBROUTINE STORE_COORDS(XSTORE,YSTORE,ZSTORE,XELM,YELM,ZELM,&
     &ISPEC,NSPEC)
              INTEGER(KIND=4), INTENT(IN) :: NSPEC
              REAL(KIND=8), INTENT(INOUT) :: XSTORE(2,2,2,NSPEC)
              REAL(KIND=8), INTENT(INOUT) :: YSTORE(2,2,2,NSPEC)
              REAL(KIND=8), INTENT(INOUT) :: ZSTORE(2,2,2,NSPEC)
              REAL(KIND=8), INTENT(IN) :: XELM(8)
              REAL(KIND=8), INTENT(IN) :: YELM(8)
              REAL(KIND=8), INTENT(IN) :: ZELM(8)
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
            END SUBROUTINE STORE_COORDS
          END INTERFACE 
        END MODULE STORE_COORDS__genmod
