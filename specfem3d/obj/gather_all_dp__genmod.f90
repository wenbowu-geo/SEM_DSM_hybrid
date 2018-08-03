        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GATHER_ALL_DP__genmod
          INTERFACE 
            SUBROUTINE GATHER_ALL_DP(SENDBUF,SENDCNT,RECVBUF,RECVCOUNT, &
     &NPROC)
              INTEGER(KIND=4) :: NPROC
              INTEGER(KIND=4) :: RECVCOUNT
              INTEGER(KIND=4) :: SENDCNT
              REAL(KIND=8) :: SENDBUF(SENDCNT)
              REAL(KIND=8) :: RECVBUF(RECVCOUNT,0:NPROC-1)
            END SUBROUTINE GATHER_ALL_DP
          END INTERFACE 
        END MODULE GATHER_ALL_DP__genmod
