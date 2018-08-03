        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:52 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE USUAL_HEX_NODES__genmod
          INTERFACE 
            SUBROUTINE USUAL_HEX_NODES(NGNOD,IADDX,IADDY,IADDZ)
              INTEGER(KIND=4), INTENT(IN) :: NGNOD
              INTEGER(KIND=4), INTENT(OUT) :: IADDX(NGNOD)
              INTEGER(KIND=4), INTENT(OUT) :: IADDY(NGNOD)
              INTEGER(KIND=4), INTENT(OUT) :: IADDZ(NGNOD)
            END SUBROUTINE USUAL_HEX_NODES
          END INTERFACE 
        END MODULE USUAL_HEX_NODES__genmod
