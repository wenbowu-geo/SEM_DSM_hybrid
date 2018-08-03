        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:50 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_JACOBIAN_2D__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_JACOBIAN_2D(MYRANK,ISPECB,XELM,YELM,ZELM,&
     &DERSHAPE2D,JACOBIAN2D,NORMAL,NGLLA,NGLLB,NSPEC2DMAX_AB)
              INTEGER(KIND=4) :: NSPEC2DMAX_AB
              INTEGER(KIND=4) :: NGLLB
              INTEGER(KIND=4) :: NGLLA
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: ISPECB
              REAL(KIND=8) :: XELM(4)
              REAL(KIND=8) :: YELM(4)
              REAL(KIND=8) :: ZELM(4)
              REAL(KIND=8) :: DERSHAPE2D(2,4,NGLLA,NGLLB)
              REAL(KIND=8) :: JACOBIAN2D(NGLLA,NGLLB,NSPEC2DMAX_AB)
              REAL(KIND=8) :: NORMAL(3,NGLLA,NGLLB,NSPEC2DMAX_AB)
            END SUBROUTINE COMPUTE_JACOBIAN_2D
          END INTERFACE 
        END MODULE COMPUTE_JACOBIAN_2D__genmod
