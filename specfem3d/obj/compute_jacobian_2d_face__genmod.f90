        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:33 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_JACOBIAN_2D_FACE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_JACOBIAN_2D_FACE(MYRANK,XELM,YELM,ZELM,  &
     &DERSHAPE2D,WGLLWGLL,JACOBIAN2DW_FACE,NORMAL_FACE,NGLLA,NGLLB,     &
     &NGNOD2D)
              INTEGER(KIND=4) :: NGNOD2D
              INTEGER(KIND=4) :: NGLLB
              INTEGER(KIND=4) :: NGLLA
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: XELM(NGNOD2D)
              REAL(KIND=8) :: YELM(NGNOD2D)
              REAL(KIND=8) :: ZELM(NGNOD2D)
              REAL(KIND=8) :: DERSHAPE2D(2,NGNOD2D,NGLLA,NGLLB)
              REAL(KIND=8) :: WGLLWGLL(NGLLA,NGLLB)
              REAL(KIND=8) :: JACOBIAN2DW_FACE(NGLLA,NGLLB)
              REAL(KIND=8) :: NORMAL_FACE(3,NGLLA,NGLLB)
            END SUBROUTINE COMPUTE_JACOBIAN_2D_FACE
          END INTERFACE 
        END MODULE COMPUTE_JACOBIAN_2D_FACE__genmod
