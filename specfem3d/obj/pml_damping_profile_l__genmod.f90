        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:20 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PML_DAMPING_PROFILE_L__genmod
          INTERFACE 
            FUNCTION PML_DAMPING_PROFILE_L(MYRANK,IGLOB,DIST,VP,DELTA)
              INTEGER(KIND=4), INTENT(IN) :: MYRANK
              INTEGER(KIND=4), INTENT(IN) :: IGLOB
              REAL(KIND=8) :: DIST
              REAL(KIND=8), INTENT(IN) :: VP
              REAL(KIND=8), INTENT(IN) :: DELTA
              REAL(KIND=8) :: PML_DAMPING_PROFILE_L
            END FUNCTION PML_DAMPING_PROFILE_L
          END INTERFACE 
        END MODULE PML_DAMPING_PROFILE_L__genmod
