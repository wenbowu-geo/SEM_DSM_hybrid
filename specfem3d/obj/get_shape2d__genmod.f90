        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:33 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_SHAPE2D__genmod
          INTERFACE 
            SUBROUTINE GET_SHAPE2D(MYRANK,SHAPE2D,DERSHAPE2D,XIGLL,YIGLL&
     &,NGLLA,NGLLB,NGNOD,NGNOD2D)
              INTEGER(KIND=4) :: NGNOD2D
              INTEGER(KIND=4) :: NGLLB
              INTEGER(KIND=4) :: NGLLA
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: SHAPE2D(NGNOD2D,NGLLA,NGLLB)
              REAL(KIND=8) :: DERSHAPE2D(2,NGNOD2D,NGLLA,NGLLB)
              REAL(KIND=8) :: XIGLL(NGLLA)
              REAL(KIND=8) :: YIGLL(NGLLB)
              INTEGER(KIND=4) :: NGNOD
            END SUBROUTINE GET_SHAPE2D
          END INTERFACE 
        END MODULE GET_SHAPE2D__genmod
