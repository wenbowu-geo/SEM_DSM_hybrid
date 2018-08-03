        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:02 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_CMT__genmod
          INTERFACE 
            SUBROUTINE GET_CMT(YR,JDA,HO,MI,SEC,TSHIFT_CMT,HDUR,LAT,LONG&
     &,DEPTH,MOMENT_TENSOR,DT,NSOURCES,MIN_TSHIFT_CMT_ORIGINAL,         &
     &USER_SOURCE_TIME_FUNCTION)
              USE SHARED_PARAMETERS, ONLY :                             &
     &          NUMBER_OF_SIMULTANEOUS_RUNS,                            &
     &          EXTERNAL_STF,                                           &
     &          NSTEP_STF,                                              &
     &          NSOURCES_STF
              INTEGER(KIND=4), INTENT(IN) :: NSOURCES
              INTEGER(KIND=4), INTENT(OUT) :: YR
              INTEGER(KIND=4), INTENT(OUT) :: JDA
              INTEGER(KIND=4), INTENT(OUT) :: HO
              INTEGER(KIND=4), INTENT(OUT) :: MI
              REAL(KIND=8), INTENT(OUT) :: SEC
              REAL(KIND=8), INTENT(OUT) :: TSHIFT_CMT(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: HDUR(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: LAT(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: LONG(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: DEPTH(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: MOMENT_TENSOR(6,NSOURCES)
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8), INTENT(OUT) :: MIN_TSHIFT_CMT_ORIGINAL
              REAL(KIND=8), INTENT(OUT) :: USER_SOURCE_TIME_FUNCTION(   &
     &NSTEP_STF,NSOURCES_STF)
            END SUBROUTINE GET_CMT
          END INTERFACE 
        END MODULE GET_CMT__genmod
