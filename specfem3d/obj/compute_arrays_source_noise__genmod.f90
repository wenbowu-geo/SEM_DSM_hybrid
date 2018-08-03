        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_ARRAYS_SOURCE_NOISE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_ARRAYS_SOURCE_NOISE(MYRANK,XI_NOISE,     &
     &ETA_NOISE,GAMMA_NOISE,NU_SINGLE,NOISE_SOURCEARRAY,XIGLL,YIGLL,    &
     &ZIGLL,NSTEP)
              INTEGER(KIND=4) :: NSTEP
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: XI_NOISE
              REAL(KIND=8) :: ETA_NOISE
              REAL(KIND=8) :: GAMMA_NOISE
              REAL(KIND=8) :: NU_SINGLE(3,3)
              REAL(KIND=8) :: NOISE_SOURCEARRAY(3,5,5,5,NSTEP)
              REAL(KIND=8) :: XIGLL(5)
              REAL(KIND=8) :: YIGLL(5)
              REAL(KIND=8) :: ZIGLL(5)
            END SUBROUTINE COMPUTE_ARRAYS_SOURCE_NOISE
          END INTERFACE 
        END MODULE COMPUTE_ARRAYS_SOURCE_NOISE__genmod
