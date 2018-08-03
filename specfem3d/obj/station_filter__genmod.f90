        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:08 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STATION_FILTER__genmod
          INTERFACE 
            SUBROUTINE STATION_FILTER(SUPPRESS_UTM_PROJECTION,          &
     &UTM_PROJECTION_ZONE,MYRANK,FILENAME,FILTERED_FILENAME,NFILTER,    &
     &LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX)
              LOGICAL(KIND=4) :: SUPPRESS_UTM_PROJECTION
              INTEGER(KIND=4) :: UTM_PROJECTION_ZONE
              INTEGER(KIND=4) :: MYRANK
              CHARACTER(*) :: FILENAME
              CHARACTER(*) :: FILTERED_FILENAME
              INTEGER(KIND=4) :: NFILTER
              REAL(KIND=8) :: LATITUDE_MIN
              REAL(KIND=8) :: LATITUDE_MAX
              REAL(KIND=8) :: LONGITUDE_MIN
              REAL(KIND=8) :: LONGITUDE_MAX
            END SUBROUTINE STATION_FILTER
          END INTERFACE 
        END MODULE STATION_FILTER__genmod
