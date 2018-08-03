        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:45:17 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SUM_KERNEL__genmod
          INTERFACE 
            SUBROUTINE SUM_KERNEL(KERNEL_NAME,KERNEL_LIST,NKER)
              CHARACTER(LEN=512) :: KERNEL_NAME
              CHARACTER(LEN=512) :: KERNEL_LIST(10000)
              INTEGER(KIND=4) :: NKER
            END SUBROUTINE SUM_KERNEL
          END INTERFACE 
        END MODULE SUM_KERNEL__genmod
