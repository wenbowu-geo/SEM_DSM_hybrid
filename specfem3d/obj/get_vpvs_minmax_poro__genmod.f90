        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:28 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_VPVS_MINMAX_PORO__genmod
          INTERFACE 
            SUBROUTINE GET_VPVS_MINMAX_PORO(VPMIN,VPMAX,VP2MIN,VP2MAX,  &
     &VSMIN,VSMAX,ISPEC,HAS_VS_ZERO,HAS_VP2_ZERO,NSPEC_AB,PHISTORE,     &
     &TORTSTORE,RHOARRAYSTORE,RHO_VPI,RHO_VPII,RHO_VSI)
              INTEGER(KIND=4) :: NSPEC_AB
              REAL(KIND=8) :: VPMIN
              REAL(KIND=8) :: VPMAX
              REAL(KIND=8) :: VP2MIN
              REAL(KIND=8) :: VP2MAX
              REAL(KIND=8) :: VSMIN
              REAL(KIND=8) :: VSMAX
              INTEGER(KIND=4) :: ISPEC
              LOGICAL(KIND=4) :: HAS_VS_ZERO
              LOGICAL(KIND=4) :: HAS_VP2_ZERO
              REAL(KIND=8) :: PHISTORE(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: TORTSTORE(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHOARRAYSTORE(2,5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHO_VPI(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHO_VPII(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHO_VSI(5,5,5,NSPEC_AB)
            END SUBROUTINE GET_VPVS_MINMAX_PORO
          END INTERFACE 
        END MODULE GET_VPVS_MINMAX_PORO__genmod
