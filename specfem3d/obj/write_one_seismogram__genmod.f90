        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:38 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_ONE_SEISMOGRAM__genmod
          INTERFACE 
            SUBROUTINE WRITE_ONE_SEISMOGRAM(ONE_SEISMOGRAM,IREC,        &
     &STATION_NAME,NETWORK_NAME,NREC,DT,T0,IT,NSTEP,SIMULATION_TYPE,    &
     &MYRANK,COMPONENT,ISTORE)
              INTEGER(KIND=4), INTENT(IN) :: NSTEP
              INTEGER(KIND=4) :: NREC
              REAL(KIND=8) :: ONE_SEISMOGRAM(3,NSTEP)
              INTEGER(KIND=4) :: IREC
              CHARACTER(LEN=32) :: STATION_NAME(NREC)
              CHARACTER(LEN=8) :: NETWORK_NAME(NREC)
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: T0
              INTEGER(KIND=4), INTENT(IN) :: IT
              INTEGER(KIND=4), INTENT(IN) :: SIMULATION_TYPE
              INTEGER(KIND=4) :: MYRANK
              CHARACTER(LEN=1) :: COMPONENT
              INTEGER(KIND=4), INTENT(IN) :: ISTORE
            END SUBROUTINE WRITE_ONE_SEISMOGRAM
          END INTERFACE 
        END MODULE WRITE_ONE_SEISMOGRAM__genmod
