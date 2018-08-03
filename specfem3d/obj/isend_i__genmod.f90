        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ISEND_I__genmod
          INTERFACE 
            SUBROUTINE ISEND_I(SENDBUF,SENDCOUNT,DEST,SENDTAG,REQ)
              INTEGER(KIND=4) :: SENDCOUNT
              INTEGER(KIND=4) :: SENDBUF(SENDCOUNT)
              INTEGER(KIND=4) :: DEST
              INTEGER(KIND=4) :: SENDTAG
              INTEGER(KIND=4) :: REQ
            END SUBROUTINE ISEND_I
          END INTERFACE 
        END MODULE ISEND_I__genmod
