        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:38 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_ADJ_SEISMOGRAMS_TO_FILE__genmod
          INTERFACE 
            SUBROUTINE WRITE_ADJ_SEISMOGRAMS_TO_FILE(MYRANK,SEISMOGRAMS,&
     &NUMBER_RECEIVER_GLOBAL,NREC_LOCAL,IT,DT,NSTEP,T0,ISTORE)
              INTEGER(KIND=4) :: NSTEP
              INTEGER(KIND=4) :: NREC_LOCAL
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: SEISMOGRAMS(3,NREC_LOCAL,NSTEP)
              INTEGER(KIND=4) :: NUMBER_RECEIVER_GLOBAL(NREC_LOCAL)
              INTEGER(KIND=4) :: IT
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: T0
              INTEGER(KIND=4) :: ISTORE
            END SUBROUTINE WRITE_ADJ_SEISMOGRAMS_TO_FILE
          END INTERFACE 
        END MODULE WRITE_ADJ_SEISMOGRAMS_TO_FILE__genmod
