        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:38 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_OUTPUT_ASCII_OR_BINARY__genmod
          INTERFACE 
            SUBROUTINE WRITE_OUTPUT_ASCII_OR_BINARY(ONE_SEISMOGRAM,NSTEP&
     &,IT,SIMULATION_TYPE,DT,T0,IORIENTATION,SISNAME,FINAL_LOCAL_PATH)
              INTEGER(KIND=4) :: NSTEP
              REAL(KIND=8) :: ONE_SEISMOGRAM(3,NSTEP)
              INTEGER(KIND=4) :: IT
              INTEGER(KIND=4) :: SIMULATION_TYPE
              REAL(KIND=8) :: DT
              REAL(KIND=8) :: T0
              INTEGER(KIND=4) :: IORIENTATION
              CHARACTER(LEN=512) :: SISNAME
              CHARACTER(LEN=512) :: FINAL_LOCAL_PATH
            END SUBROUTINE WRITE_OUTPUT_ASCII_OR_BINARY
          END INTERFACE 
        END MODULE WRITE_OUTPUT_ASCII_OR_BINARY__genmod
