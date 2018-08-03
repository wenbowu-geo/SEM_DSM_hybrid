        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:34 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RECOMPUTE_JACOBIAN_27__genmod
          INTERFACE 
            SUBROUTINE RECOMPUTE_JACOBIAN_27(NGNOD,NDIM,XI,ETA,GAMMA,   &
     &SHAPE3D,DERSHAPE3D)
              INTEGER(KIND=4) :: NDIM
              INTEGER(KIND=4) :: NGNOD
              REAL(KIND=8) :: XI
              REAL(KIND=8) :: ETA
              REAL(KIND=8) :: GAMMA
              REAL(KIND=8) :: SHAPE3D(NGNOD)
              REAL(KIND=8) :: DERSHAPE3D(NDIM,NGNOD)
            END SUBROUTINE RECOMPUTE_JACOBIAN_27
          END INTERFACE 
        END MODULE RECOMPUTE_JACOBIAN_27__genmod
