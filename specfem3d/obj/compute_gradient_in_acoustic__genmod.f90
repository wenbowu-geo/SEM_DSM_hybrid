        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:58 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_GRADIENT_IN_ACOUSTIC__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_GRADIENT_IN_ACOUSTIC(ISPEC,NSPEC_AB,     &
     &NGLOB_AB,SCALAR_FIELD,VECTOR_FIELD_ELEMENT,HPRIME_XX,HPRIME_YY,   &
     &HPRIME_ZZ,XIX,XIY,XIZ,ETAX,ETAY,ETAZ,GAMMAX,GAMMAY,GAMMAZ,IBOOL,  &
     &RHOSTORE,GRAVITY)
              INTEGER(KIND=4), INTENT(IN) :: NGLOB_AB
              INTEGER(KIND=4), INTENT(IN) :: NSPEC_AB
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
              REAL(KIND=8), INTENT(IN) :: SCALAR_FIELD(NGLOB_AB)
              REAL(KIND=8), INTENT(OUT) :: VECTOR_FIELD_ELEMENT(3,5,5,5)
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
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHOSTORE(5,5,5,NSPEC_AB)
              LOGICAL(KIND=4) :: GRAVITY
            END SUBROUTINE COMPUTE_GRADIENT_IN_ACOUSTIC
          END INTERFACE 
        END MODULE COMPUTE_GRADIENT_IN_ACOUSTIC__genmod
