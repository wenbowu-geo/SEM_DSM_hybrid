        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:37 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_ACOUSTIC_LDDRK__genmod
          INTERFACE 
            SUBROUTINE UPDATE_ACOUSTIC_LDDRK(NGLOB,NGLOB_LDDRK,         &
     &POTENTIAL_ACOUSTIC,POTENTIAL_DOT_ACOUSTIC,                        &
     &POTENTIAL_DOT_DOT_ACOUSTIC,POTENTIAL_ACOUSTIC_LDDRK,              &
     &POTENTIAL_DOT_ACOUSTIC_LDDRK,DELTAT,ALPHA,BETA)
              INTEGER(KIND=4), INTENT(IN) :: NGLOB_LDDRK
              INTEGER(KIND=4), INTENT(IN) :: NGLOB
              REAL(KIND=8), INTENT(INOUT) :: POTENTIAL_ACOUSTIC(NGLOB)
              REAL(KIND=8), INTENT(INOUT) :: POTENTIAL_DOT_ACOUSTIC(    &
     &NGLOB)
              REAL(KIND=8), INTENT(INOUT) :: POTENTIAL_DOT_DOT_ACOUSTIC(&
     &NGLOB)
              REAL(KIND=8), INTENT(INOUT) :: POTENTIAL_ACOUSTIC_LDDRK(  &
     &NGLOB_LDDRK)
              REAL(KIND=8), INTENT(INOUT) ::                            &
     &POTENTIAL_DOT_ACOUSTIC_LDDRK(NGLOB_LDDRK)
              REAL(KIND=8), INTENT(IN) :: DELTAT
              REAL(KIND=8), INTENT(IN) :: ALPHA
              REAL(KIND=8), INTENT(IN) :: BETA
            END SUBROUTINE UPDATE_ACOUSTIC_LDDRK
          END INTERFACE 
        END MODULE UPDATE_ACOUSTIC_LDDRK__genmod
