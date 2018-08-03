        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:45 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CVD_WRITE_GLL_POINTS__genmod
          INTERFACE 
            SUBROUTINE CVD_WRITE_GLL_POINTS(NSPEC_AB,NGLOB_AB,IBOOL,    &
     &XSTORE,YSTORE,ZSTORE,DAT,IT,NPP,NUMPOIN,NP)
              INTEGER(KIND=4), INTENT(IN) :: NGLOB_AB
              INTEGER(KIND=4), INTENT(IN) :: NSPEC_AB
              INTEGER(KIND=4), INTENT(IN) :: IBOOL(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: XSTORE(NGLOB_AB)
              REAL(KIND=8) :: YSTORE(NGLOB_AB)
              REAL(KIND=8) :: ZSTORE(NGLOB_AB)
              REAL(KIND=4), INTENT(IN) :: DAT(5,5,5,NSPEC_AB)
              INTEGER(KIND=4) :: IT
              INTEGER(KIND=4) :: NPP
              INTEGER(KIND=4) :: NUMPOIN
              INTEGER(KIND=4) :: NP
            END SUBROUTINE CVD_WRITE_GLL_POINTS
          END INTERFACE 
        END MODULE CVD_WRITE_GLL_POINTS__genmod
