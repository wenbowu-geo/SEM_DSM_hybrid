        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ALL_GATHER_ALL_I__genmod
          INTERFACE 
            SUBROUTINE ALL_GATHER_ALL_I(SENDBUF,SENDCNT,RECVBUF,        &
     &RECVCOUNT,NPROC)
              INTEGER(KIND=4) :: NPROC
              INTEGER(KIND=4) :: RECVCOUNT
              INTEGER(KIND=4) :: SENDCNT
              INTEGER(KIND=4) :: SENDBUF(SENDCNT)
              INTEGER(KIND=4) :: RECVBUF(RECVCOUNT,0:NPROC-1)
            END SUBROUTINE ALL_GATHER_ALL_I
          END INTERFACE 
        END MODULE ALL_GATHER_ALL_I__genmod
