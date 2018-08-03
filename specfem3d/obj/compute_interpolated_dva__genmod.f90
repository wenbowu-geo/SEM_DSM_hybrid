        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:58 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_INTERPOLATED_DVA__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_INTERPOLATED_DVA(DISPL,VELOC,ACCEL,      &
     &NGLOB_AB,ISPEC,NSPEC_AB,IBOOL,XI_R,ETA_R,GAMMA_R,HXIR,HETAR,      &
     &HGAMMAR,DXD,DYD,DZD,VXD,VYD,VZD,AXD,AYD,AZD)
              INTEGER(KIND=4), INTENT(IN) :: NSPEC_AB
              INTEGER(KIND=4), INTENT(IN) :: NGLOB_AB
              REAL(KIND=8), INTENT(IN) :: DISPL(3,NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: VELOC(3,NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: ACCEL(3,NGLOB_AB)
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
              INTEGER(KIND=4), INTENT(IN) :: IBOOL(5,5,5,NSPEC_AB)
              REAL(KIND=8), INTENT(IN) :: XI_R
              REAL(KIND=8), INTENT(IN) :: ETA_R
              REAL(KIND=8), INTENT(IN) :: GAMMA_R
              REAL(KIND=8), INTENT(IN) :: HXIR(5)
              REAL(KIND=8), INTENT(IN) :: HETAR(5)
              REAL(KIND=8), INTENT(IN) :: HGAMMAR(5)
              REAL(KIND=8), INTENT(OUT) :: DXD
              REAL(KIND=8), INTENT(OUT) :: DYD
              REAL(KIND=8), INTENT(OUT) :: DZD
              REAL(KIND=8), INTENT(OUT) :: VXD
              REAL(KIND=8), INTENT(OUT) :: VYD
              REAL(KIND=8), INTENT(OUT) :: VZD
              REAL(KIND=8), INTENT(OUT) :: AXD
              REAL(KIND=8), INTENT(OUT) :: AYD
              REAL(KIND=8), INTENT(OUT) :: AZD
            END SUBROUTINE COMPUTE_INTERPOLATED_DVA
          END INTERFACE 
        END MODULE COMPUTE_INTERPOLATED_DVA__genmod
