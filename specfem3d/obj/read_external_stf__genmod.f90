        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:24 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_EXTERNAL_STF__genmod
          INTERFACE 
            SUBROUTINE READ_EXTERNAL_STF(ISOURCE,                       &
     &USER_SOURCE_TIME_FUNCTION,EXTERNAL_STF_FILENAME)
              USE SHARED_PARAMETERS, ONLY :                             &
     &          NSTEP,                                                  &
     &          DT,                                                     &
     &          NSTEP_STF,                                              &
     &          NSOURCES_STF
              INTEGER(KIND=4), INTENT(IN) :: ISOURCE
              REAL(KIND=8), INTENT(INOUT) :: USER_SOURCE_TIME_FUNCTION( &
     &NSTEP_STF,NSOURCES_STF)
              CHARACTER(LEN=512), INTENT(IN) :: EXTERNAL_STF_FILENAME
            END SUBROUTINE READ_EXTERNAL_STF
          END INTERFACE 
        END MODULE READ_EXTERNAL_STF__genmod
