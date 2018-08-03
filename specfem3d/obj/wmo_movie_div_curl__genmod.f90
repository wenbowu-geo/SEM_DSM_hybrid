        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:37 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WMO_MOVIE_DIV_CURL__genmod
          INTERFACE 
            SUBROUTINE WMO_MOVIE_DIV_CURL(NSPEC_AB,NGLOB_AB,VELOC,      &
     &DIV_GLOB,VALENCE,DIV,CURL_X,CURL_Y,CURL_Z,VELOCITY_X,VELOCITY_Y,  &
     &VELOCITY_Z,IBOOL,ISPEC_IS,HPRIME_XX,HPRIME_YY,HPRIME_ZZ,XIX,XIY,  &
     &XIZ,ETAX,ETAY,ETAZ,GAMMAX,GAMMAY,GAMMAZ)
              INTEGER(KIND=4) :: NGLOB_AB
              INTEGER(KIND=4) :: NSPEC_AB
              REAL(KIND=8), INTENT(IN) :: VELOC(3,NGLOB_AB)
              REAL(KIND=8) :: DIV_GLOB(NGLOB_AB)
              INTEGER(KIND=4) :: VALENCE(NGLOB_AB)
              REAL(KIND=8) :: DIV(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: CURL_X(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: CURL_Y(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: CURL_Z(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: VELOCITY_X(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: VELOCITY_Y(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: VELOCITY_Z(5,5,5,NSPEC_AB)
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC_AB)
              LOGICAL(KIND=4) :: ISPEC_IS(NSPEC_AB)
              REAL(KIND=8) :: HPRIME_XX(5,5)
              REAL(KIND=8) :: HPRIME_YY(5,5)
              REAL(KIND=8) :: HPRIME_ZZ(5,5)
              REAL(KIND=8) :: XIX(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: XIY(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: XIZ(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: ETAX(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: ETAY(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: ETAZ(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: GAMMAX(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: GAMMAY(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: GAMMAZ(5,5,5,NSPEC_AB)
            END SUBROUTINE WMO_MOVIE_DIV_CURL
          END INTERFACE 
        END MODULE WMO_MOVIE_DIV_CURL__genmod
