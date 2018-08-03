        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SUM_ALL_1DARRAY_DP__genmod
          INTERFACE 
            SUBROUTINE SUM_ALL_1DARRAY_DP(SENDBUF,RECVBUF,NX)
              INTEGER(KIND=4) :: NX
              REAL(KIND=8) :: SENDBUF(NX)
              REAL(KIND=8) :: RECVBUF(NX)
            END SUBROUTINE SUM_ALL_1DARRAY_DP
          END INTERFACE 
        END MODULE SUM_ALL_1DARRAY_DP__genmod
