        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:34 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RECOMPUTE_JACOBIAN__genmod
          INTERFACE 
            SUBROUTINE RECOMPUTE_JACOBIAN(XELM,YELM,ZELM,XI,ETA,GAMMA,X,&
     &Y,Z,XIX,XIY,XIZ,ETAX,ETAY,ETAZ,GAMMAX,GAMMAY,GAMMAZ,NGNOD)
              INTEGER(KIND=4) :: NGNOD
              REAL(KIND=8) :: XELM(NGNOD)
              REAL(KIND=8) :: YELM(NGNOD)
              REAL(KIND=8) :: ZELM(NGNOD)
              REAL(KIND=8) :: XI
              REAL(KIND=8) :: ETA
              REAL(KIND=8) :: GAMMA
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y
              REAL(KIND=8) :: Z
              REAL(KIND=8) :: XIX
              REAL(KIND=8) :: XIY
              REAL(KIND=8) :: XIZ
              REAL(KIND=8) :: ETAX
              REAL(KIND=8) :: ETAY
              REAL(KIND=8) :: ETAZ
              REAL(KIND=8) :: GAMMAX
              REAL(KIND=8) :: GAMMAY
              REAL(KIND=8) :: GAMMAZ
            END SUBROUTINE RECOMPUTE_JACOBIAN
          END INTERFACE 
        END MODULE RECOMPUTE_JACOBIAN__genmod
