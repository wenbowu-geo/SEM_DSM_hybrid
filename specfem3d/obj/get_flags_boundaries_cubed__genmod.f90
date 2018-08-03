        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:48 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_FLAGS_BOUNDARIES_CUBED__genmod
          INTERFACE 
            SUBROUTINE GET_FLAGS_BOUNDARIES_CUBED(NSPEC,IPROC_XI,       &
     &IPROC_ETA,ISPEC,IDOUBLING,XELM,YELM,ZELM,IBOUN,IMPICUT_XI,        &
     &IMPICUT_ETA,NPROC_XI,NPROC_ETA,ANGULAR_WIDTH_XI_IN_DEGREES,       &
     &ANGULAR_WIDTH_ETA_IN_DEGREES,Z_DEPTH_BLOCK,                       &
     &IFLAG_ONE_LAYER_TOPOGRAPHY)
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: IPROC_XI
              INTEGER(KIND=4) :: IPROC_ETA
              INTEGER(KIND=4) :: ISPEC
              INTEGER(KIND=4) :: IDOUBLING
              REAL(KIND=8) :: XELM(8)
              REAL(KIND=8) :: YELM(8)
              REAL(KIND=8) :: ZELM(8)
              LOGICAL(KIND=4) :: IBOUN(6,NSPEC)
              LOGICAL(KIND=4) :: IMPICUT_XI(2,NSPEC)
              LOGICAL(KIND=4) :: IMPICUT_ETA(2,NSPEC)
              INTEGER(KIND=4) :: NPROC_XI
              INTEGER(KIND=4) :: NPROC_ETA
              REAL(KIND=8) :: ANGULAR_WIDTH_XI_IN_DEGREES
              REAL(KIND=8) :: ANGULAR_WIDTH_ETA_IN_DEGREES
              REAL(KIND=8) :: Z_DEPTH_BLOCK
              INTEGER(KIND=4) :: IFLAG_ONE_LAYER_TOPOGRAPHY
            END SUBROUTINE GET_FLAGS_BOUNDARIES_CUBED
          END INTERFACE 
        END MODULE GET_FLAGS_BOUNDARIES_CUBED__genmod
