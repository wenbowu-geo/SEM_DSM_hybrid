        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RECV_I__genmod
          INTERFACE 
            SUBROUTINE RECV_I(RECVBUF,RECVCOUNT,DEST,RECVTAG)
              INTEGER(KIND=4) :: RECVCOUNT
              INTEGER(KIND=4) :: RECVBUF(RECVCOUNT)
              INTEGER(KIND=4) :: DEST
              INTEGER(KIND=4) :: RECVTAG
            END SUBROUTINE RECV_I
          END INTERFACE 
        END MODULE RECV_I__genmod
