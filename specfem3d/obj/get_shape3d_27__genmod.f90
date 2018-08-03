        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:51 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_SHAPE3D_27__genmod
          INTERFACE 
            SUBROUTINE GET_SHAPE3D_27(NGNOD,SHAPE3D,DERSHAPE3D,XI,ETA,  &
     &GAMMA,I,J,K)
              INTEGER(KIND=4) :: NGNOD
              REAL(KIND=8) :: SHAPE3D(NGNOD,5,5,5)
              REAL(KIND=8) :: DERSHAPE3D(3,NGNOD,5,5,5)
              REAL(KIND=8) :: XI
              REAL(KIND=8) :: ETA
              REAL(KIND=8) :: GAMMA
              INTEGER(KIND=4) :: I
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: K
            END SUBROUTINE GET_SHAPE3D_27
          END INTERFACE 
        END MODULE GET_SHAPE3D_27__genmod
