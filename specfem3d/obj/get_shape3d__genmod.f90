        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:51 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_SHAPE3D__genmod
          INTERFACE 
            SUBROUTINE GET_SHAPE3D(MYRANK,SHAPE3D,DERSHAPE3D,XIGLL,YIGLL&
     &,ZIGLL,NGNOD)
              INTEGER(KIND=4) :: NGNOD
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: SHAPE3D(NGNOD,5,5,5)
              REAL(KIND=8) :: DERSHAPE3D(3,NGNOD,5,5,5)
              REAL(KIND=8) :: XIGLL(5)
              REAL(KIND=8) :: YIGLL(5)
              REAL(KIND=8) :: ZIGLL(5)
            END SUBROUTINE GET_SHAPE3D
          END INTERFACE 
        END MODULE GET_SHAPE3D__genmod
