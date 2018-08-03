        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MAX_ALLREDUCE_I__genmod
          INTERFACE 
            SUBROUTINE MAX_ALLREDUCE_I(BUFFER,COUNTVAL)
              INTEGER(KIND=4) :: COUNTVAL
              INTEGER(KIND=4), INTENT(INOUT) :: BUFFER(COUNTVAL)
            END SUBROUTINE MAX_ALLREDUCE_I
          END INTERFACE 
        END MODULE MAX_ALLREDUCE_I__genmod
