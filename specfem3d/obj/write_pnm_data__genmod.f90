        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:01 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_PNM_DATA__genmod
          INTERFACE 
            SUBROUTINE WRITE_PNM_DATA(COLOR_IMAGE_2D_DATA,              &
     &IGLOB_IMAGE_COLOR_2D,NX,NY,IT,CUTSNAPS,IMAGE_COLOR_VP_DISPLAY)
              INTEGER(KIND=4), INTENT(IN) :: NY
              INTEGER(KIND=4), INTENT(IN) :: NX
              REAL(KIND=8) :: COLOR_IMAGE_2D_DATA(NX,NY)
              INTEGER(KIND=4) :: IGLOB_IMAGE_COLOR_2D(NX,NY)
              INTEGER(KIND=4), INTENT(IN) :: IT
              REAL(KIND=8) :: CUTSNAPS
              REAL(KIND=8) :: IMAGE_COLOR_VP_DISPLAY(NX,NY)
            END SUBROUTINE WRITE_PNM_DATA
          END INTERFACE 
        END MODULE WRITE_PNM_DATA__genmod
