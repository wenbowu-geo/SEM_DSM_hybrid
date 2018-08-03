        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:38 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_SEISMOGRAMS_TO_FILE__genmod
          INTERFACE 
            SUBROUTINE WRITE_SEISMOGRAMS_TO_FILE(SEISMOGRAMS,ISTORE)
              USE SPECFEM_PAR, ONLY :                                   &
     &          MYRANK,                                                 &
     &          NUMBER_RECEIVER_GLOBAL,                                 &
     &          STATION_NAME,                                           &
     &          NETWORK_NAME,                                           &
     &          NREC,                                                   &
     &          NREC_LOCAL,                                             &
     &          ISLICE_SELECTED_REC,                                    &
     &          IT,                                                     &
     &          DT,                                                     &
     &          NSTEP,                                                  &
     &          T0,                                                     &
     &          SIMULATION_TYPE,                                        &
     &          WRITE_SEISMOGRAMS_BY_MASTER,                            &
     &          SAVE_ALL_SEISMOS_IN_ONE_FILE,                           &
     &          USE_BINARY_FOR_SEISMOGRAMS
              REAL(KIND=8) :: SEISMOGRAMS(3,NREC_LOCAL,NSTEP)
              INTEGER(KIND=4) :: ISTORE
            END SUBROUTINE WRITE_SEISMOGRAMS_TO_FILE
          END INTERFACE 
        END MODULE WRITE_SEISMOGRAMS_TO_FILE__genmod
