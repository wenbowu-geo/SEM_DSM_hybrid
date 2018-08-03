        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:50 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMBINE_SEM_ARRAY__genmod
          INTERFACE 
            SUBROUTINE COMBINE_SEM_ARRAY(KERNEL_NAME,KERNEL_PATHS,      &
     &OUTPUT_DIR,NPATH)
              CHARACTER(LEN=512) :: KERNEL_NAME
              CHARACTER(LEN=512) :: KERNEL_PATHS(65535)
              CHARACTER(LEN=512) :: OUTPUT_DIR
              INTEGER(KIND=4) :: NPATH
            END SUBROUTINE COMBINE_SEM_ARRAY
          END INTERFACE 
        END MODULE COMBINE_SEM_ARRAY__genmod
