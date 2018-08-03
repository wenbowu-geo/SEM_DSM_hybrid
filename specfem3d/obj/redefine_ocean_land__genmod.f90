        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:44 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE REDEFINE_OCEAN_LAND__genmod
          INTERFACE 
            SUBROUTINE REDEFINE_OCEAN_LAND(XELM,YELM,ZELM,              &
     &IMATERIAL_NUMBER,ISUBREGION,NMESHREGIONS,IR,IR1,IR2,IY,IY1,IY2,IX,&
     &IX1,IX2,IPROC_XI,NPROC_XI,IPROC_ETA,NPROC_ETA,NX_NOTOPO,NY_NOTOPO)
              REAL(KIND=8) :: XELM(8)
              REAL(KIND=8) :: YELM(8)
              REAL(KIND=8) :: ZELM(8)
              INTEGER(KIND=4) :: IMATERIAL_NUMBER
              INTEGER(KIND=4) :: ISUBREGION
              INTEGER(KIND=4) :: NMESHREGIONS
              INTEGER(KIND=4) :: IR
              INTEGER(KIND=4) :: IR1
              INTEGER(KIND=4) :: IR2
              INTEGER(KIND=4) :: IY
              INTEGER(KIND=4) :: IY1
              INTEGER(KIND=4) :: IY2
              INTEGER(KIND=4) :: IX
              INTEGER(KIND=4) :: IX1
              INTEGER(KIND=4) :: IX2
              INTEGER(KIND=4) :: IPROC_XI
              INTEGER(KIND=4) :: NPROC_XI
              INTEGER(KIND=4) :: IPROC_ETA
              INTEGER(KIND=4) :: NPROC_ETA
              INTEGER(KIND=4) :: NX_NOTOPO
              INTEGER(KIND=4) :: NY_NOTOPO
            END SUBROUTINE REDEFINE_OCEAN_LAND
          END INTERFACE 
        END MODULE REDEFINE_OCEAN_LAND__genmod
