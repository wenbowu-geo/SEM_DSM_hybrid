        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MAKE_GRAVITY__genmod
          INTERFACE 
            SUBROUTINE MAKE_GRAVITY(NSPL,RSPL,GSPL,GSPL2,ROCEAN,        &
     &RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,&
     &RCMB,RICB)
              INTEGER(KIND=4) :: NSPL
              REAL(KIND=8) :: RSPL(640)
              REAL(KIND=8) :: GSPL(640)
              REAL(KIND=8) :: GSPL2(640)
              REAL(KIND=8) :: ROCEAN
              REAL(KIND=8) :: RMIDDLE_CRUST
              REAL(KIND=8) :: RMOHO
              REAL(KIND=8) :: R80
              REAL(KIND=8) :: R220
              REAL(KIND=8) :: R400
              REAL(KIND=8) :: R600
              REAL(KIND=8) :: R670
              REAL(KIND=8) :: R771
              REAL(KIND=8) :: RTOPDDOUBLEPRIME
              REAL(KIND=8) :: RCMB
              REAL(KIND=8) :: RICB
            END SUBROUTINE MAKE_GRAVITY
          END INTERFACE 
        END MODULE MAKE_GRAVITY__genmod
