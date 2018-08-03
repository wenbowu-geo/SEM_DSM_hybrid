        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:03 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_FORCE__genmod
          INTERFACE 
            SUBROUTINE GET_FORCE(TSHIFT_FORCE,HDUR,LAT,LONG,DEPTH,      &
     &NSOURCES,MIN_TSHIFT_FORCE_ORIGINAL,FACTOR_FORCE_SOURCE,           &
     &COMP_DIR_VECT_SOURCE_E,COMP_DIR_VECT_SOURCE_N,                    &
     &COMP_DIR_VECT_SOURCE_Z_UP,USER_SOURCE_TIME_FUNCTION)
              USE SHARED_PARAMETERS, ONLY :                             &
     &          NUMBER_OF_SIMULTANEOUS_RUNS,                            &
     &          EXTERNAL_STF,                                           &
     &          NSTEP_STF,                                              &
     &          NSOURCES_STF
              INTEGER(KIND=4), INTENT(IN) :: NSOURCES
              REAL(KIND=8), INTENT(OUT) :: TSHIFT_FORCE(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: HDUR(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: LAT(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: LONG(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: DEPTH(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: MIN_TSHIFT_FORCE_ORIGINAL
              REAL(KIND=8), INTENT(OUT) :: FACTOR_FORCE_SOURCE(NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: COMP_DIR_VECT_SOURCE_E(      &
     &NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: COMP_DIR_VECT_SOURCE_N(      &
     &NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: COMP_DIR_VECT_SOURCE_Z_UP(   &
     &NSOURCES)
              REAL(KIND=8), INTENT(OUT) :: USER_SOURCE_TIME_FUNCTION(   &
     &NSTEP_STF,NSOURCES_STF)
            END SUBROUTINE GET_FORCE
          END INTERFACE 
        END MODULE GET_FORCE__genmod
