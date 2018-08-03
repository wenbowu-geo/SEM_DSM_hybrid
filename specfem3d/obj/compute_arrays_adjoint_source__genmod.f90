        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:42 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_ARRAYS_ADJOINT_SOURCE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_ARRAYS_ADJOINT_SOURCE(MYRANK,            &
     &ADJ_SOURCE_FILE,XI_RECEIVER,ETA_RECEIVER,GAMMA_RECEIVER,          &
     &ADJ_SOURCEARRAY,XIGLL,YIGLL,ZIGLL,IT_SUB_ADJ,NSTEP,               &
     &NTSTEP_BETWEEN_READ_ADJSRC)
              INTEGER(KIND=4) :: NTSTEP_BETWEEN_READ_ADJSRC
              INTEGER(KIND=4) :: MYRANK
              CHARACTER(*) :: ADJ_SOURCE_FILE
              REAL(KIND=8) :: XI_RECEIVER
              REAL(KIND=8) :: ETA_RECEIVER
              REAL(KIND=8) :: GAMMA_RECEIVER
              REAL(KIND=8) :: ADJ_SOURCEARRAY(NTSTEP_BETWEEN_READ_ADJSRC&
     &,3,5,5,5)
              REAL(KIND=8) :: XIGLL(5)
              REAL(KIND=8) :: YIGLL(5)
              REAL(KIND=8) :: ZIGLL(5)
              INTEGER(KIND=4) :: IT_SUB_ADJ
              INTEGER(KIND=4) :: NSTEP
            END SUBROUTINE COMPUTE_ARRAYS_ADJOINT_SOURCE
          END INTERFACE 
        END MODULE COMPUTE_ARRAYS_ADJOINT_SOURCE__genmod
