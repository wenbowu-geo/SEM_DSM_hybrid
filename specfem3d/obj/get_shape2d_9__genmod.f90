        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:33 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_SHAPE2D_9__genmod
          INTERFACE 
            SUBROUTINE GET_SHAPE2D_9(NGNOD2D,SHAPE2D,DERSHAPE2D,XIGLL,  &
     &YIGLL,NGLLA,NGLLB)
              INTEGER(KIND=4) :: NGLLB
              INTEGER(KIND=4) :: NGLLA
              INTEGER(KIND=4) :: NGNOD2D
              REAL(KIND=8) :: SHAPE2D(NGNOD2D,NGLLA,NGLLB)
              REAL(KIND=8) :: DERSHAPE2D(2,NGNOD2D,NGLLA,NGLLB)
              REAL(KIND=8) :: XIGLL(NGLLA)
              REAL(KIND=8) :: YIGLL(NGLLB)
            END SUBROUTINE GET_SHAPE2D_9
          END INTERFACE 
        END MODULE GET_SHAPE2D_9__genmod
