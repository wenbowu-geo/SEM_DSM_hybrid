        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:42 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_ARRAYS_SOURCE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_ARRAYS_SOURCE(ISPEC_SELECTED_SOURCE,     &
     &SOURCEARRAY,MXX,MYY,MZZ,MXY,MXZ,MYZ,XIX,XIY,XIZ,ETAX,ETAY,ETAZ,   &
     &GAMMAX,GAMMAY,GAMMAZ,HXIS,HPXIS,HETAS,HPETAS,HGAMMAS,HPGAMMAS,    &
     &NSPEC)
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: ISPEC_SELECTED_SOURCE
              REAL(KIND=8) :: SOURCEARRAY(3,5,5,5)
              REAL(KIND=8) :: MXX
              REAL(KIND=8) :: MYY
              REAL(KIND=8) :: MZZ
              REAL(KIND=8) :: MXY
              REAL(KIND=8) :: MXZ
              REAL(KIND=8) :: MYZ
              REAL(KIND=8) :: XIX(5,5,5,NSPEC)
              REAL(KIND=8) :: XIY(5,5,5,NSPEC)
              REAL(KIND=8) :: XIZ(5,5,5,NSPEC)
              REAL(KIND=8) :: ETAX(5,5,5,NSPEC)
              REAL(KIND=8) :: ETAY(5,5,5,NSPEC)
              REAL(KIND=8) :: ETAZ(5,5,5,NSPEC)
              REAL(KIND=8) :: GAMMAX(5,5,5,NSPEC)
              REAL(KIND=8) :: GAMMAY(5,5,5,NSPEC)
              REAL(KIND=8) :: GAMMAZ(5,5,5,NSPEC)
              REAL(KIND=8) :: HXIS(5)
              REAL(KIND=8) :: HPXIS(5)
              REAL(KIND=8) :: HETAS(5)
              REAL(KIND=8) :: HPETAS(5)
              REAL(KIND=8) :: HGAMMAS(5)
              REAL(KIND=8) :: HPGAMMAS(5)
            END SUBROUTINE COMPUTE_ARRAYS_SOURCE
          END INTERFACE 
        END MODULE COMPUTE_ARRAYS_SOURCE__genmod
