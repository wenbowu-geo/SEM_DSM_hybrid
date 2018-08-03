        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PREM_DENSITY__genmod
          INTERFACE 
            SUBROUTINE PREM_DENSITY(X,RHO,RICB,RCMB,RTOPDDOUBLEPRIME,   &
     &R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
              REAL(KIND=8) :: X
              REAL(KIND=8) :: RHO
              REAL(KIND=8) :: RICB
              REAL(KIND=8) :: RCMB
              REAL(KIND=8) :: RTOPDDOUBLEPRIME
              REAL(KIND=8) :: R600
              REAL(KIND=8) :: R670
              REAL(KIND=8) :: R220
              REAL(KIND=8) :: R771
              REAL(KIND=8) :: R400
              REAL(KIND=8) :: R80
              REAL(KIND=8) :: RMOHO
              REAL(KIND=8) :: RMIDDLE_CRUST
              REAL(KIND=8) :: ROCEAN
            END SUBROUTINE PREM_DENSITY
          END INTERFACE 
        END MODULE PREM_DENSITY__genmod
