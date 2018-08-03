        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:32 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FMINSEARCH__genmod
          INTERFACE 
            SUBROUTINE FMINSEARCH(FUNK,X,N,ITERCOUNT,TOLF,PRNT,ERR)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: FUNK
              EXTERNAL FUNK
              REAL(KIND=8) :: X(N)
              INTEGER(KIND=4) :: ITERCOUNT
              REAL(KIND=8) :: TOLF
              INTEGER(KIND=4) :: PRNT
              INTEGER(KIND=4) :: ERR
            END SUBROUTINE FMINSEARCH
          END INTERFACE 
        END MODULE FMINSEARCH__genmod
