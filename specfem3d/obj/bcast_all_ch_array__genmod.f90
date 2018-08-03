        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BCAST_ALL_CH_ARRAY__genmod
          INTERFACE 
            SUBROUTINE BCAST_ALL_CH_ARRAY(BUFFER,COUNTVAL)
              INTEGER(KIND=4) :: COUNTVAL
              CHARACTER(LEN=512) :: BUFFER(512*COUNTVAL)
            END SUBROUTINE BCAST_ALL_CH_ARRAY
          END INTERFACE 
        END MODULE BCAST_ALL_CH_ARRAY__genmod
