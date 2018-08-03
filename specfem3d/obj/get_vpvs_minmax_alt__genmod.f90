        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:28 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_VPVS_MINMAX_ALT__genmod
          INTERFACE 
            SUBROUTINE GET_VPVS_MINMAX_ALT(VPMIN,VPMAX,VSMIN,VSMAX,ISPEC&
     &,HAS_VS_ZERO,NSPEC_AB,KAPPASTORE,MUSTORE,RHO_VP,RHO_VS)
              INTEGER(KIND=4) :: NSPEC_AB
              REAL(KIND=8) :: VPMIN
              REAL(KIND=8) :: VPMAX
              REAL(KIND=8) :: VSMIN
              REAL(KIND=8) :: VSMAX
              INTEGER(KIND=4) :: ISPEC
              LOGICAL(KIND=4) :: HAS_VS_ZERO
              REAL(KIND=8) :: KAPPASTORE(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: MUSTORE(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHO_VP(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: RHO_VS(5,5,5,NSPEC_AB)
            END SUBROUTINE GET_VPVS_MINMAX_ALT
          END INTERFACE 
        END MODULE GET_VPVS_MINMAX_ALT__genmod
