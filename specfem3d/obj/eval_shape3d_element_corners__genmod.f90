        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:51 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EVAL_SHAPE3D_ELEMENT_CORNERS__genmod
          INTERFACE 
            SUBROUTINE EVAL_SHAPE3D_ELEMENT_CORNERS(XELM,YELM,ZELM,ISPEC&
     &,IBOOL,XSTORE,YSTORE,ZSTORE,NSPEC_AB,NGLOB_AB)
              INTEGER(KIND=4) :: NGLOB_AB
              INTEGER(KIND=4) :: NSPEC_AB
              REAL(KIND=8), INTENT(OUT) :: XELM(8)
              REAL(KIND=8), INTENT(OUT) :: YELM(8)
              REAL(KIND=8), INTENT(OUT) :: ZELM(8)
              INTEGER(KIND=4) :: ISPEC
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: XSTORE(NGLOB_AB)
              REAL(KIND=8) :: YSTORE(NGLOB_AB)
              REAL(KIND=8) :: ZSTORE(NGLOB_AB)
            END SUBROUTINE EVAL_SHAPE3D_ELEMENT_CORNERS
          END INTERFACE 
        END MODULE EVAL_SHAPE3D_ELEMENT_CORNERS__genmod
