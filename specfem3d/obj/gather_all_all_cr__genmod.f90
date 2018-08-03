        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GATHER_ALL_ALL_CR__genmod
          INTERFACE 
            SUBROUTINE GATHER_ALL_ALL_CR(SENDBUF,RECVBUF,COUNTS,NPROC)
              INTEGER(KIND=4) :: NPROC
              INTEGER(KIND=4) :: COUNTS
              REAL(KIND=8) :: SENDBUF(COUNTS)
              REAL(KIND=8) :: RECVBUF(COUNTS,0:NPROC-1)
            END SUBROUTINE GATHER_ALL_ALL_CR
          END INTERFACE 
        END MODULE GATHER_ALL_ALL_CR__genmod
