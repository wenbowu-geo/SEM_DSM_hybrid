        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PML_COMPUTE_ACCEL_CONTRIBUTION_ACOUSTIC__genmod
          INTERFACE 
            SUBROUTINE PML_COMPUTE_ACCEL_CONTRIBUTION_ACOUSTIC(ISPEC,   &
     &ISPEC_CPML,POTENTIAL_ACOUSTIC,POTENTIAL_DOT_ACOUSTIC,             &
     &RMEMORY_POTENTIAL_ACOUSTIC)
              USE SPECFEM_PAR, ONLY :                                   &
     &          NGLOB_AB,                                               &
     &          DELTAT,                                                 &
     &          WGLL_CUBE,                                              &
     &          JACOBIAN,                                               &
     &          IBOOL,                                                  &
     &          KAPPASTORE
              USE PML_PAR, ONLY :                                       &
     &          CPML_REGIONS,                                           &
     &          NSPEC_CPML,                                             &
     &          D_STORE_X,                                              &
     &          D_STORE_Y,                                              &
     &          D_STORE_Z,                                              &
     &          K_STORE_X,                                              &
     &          K_STORE_Y,                                              &
     &          K_STORE_Z,                                              &
     &          ALPHA_STORE_X,                                          &
     &          ALPHA_STORE_Y,                                          &
     &          ALPHA_STORE_Z,                                          &
     &          POTENTIAL_DOT_DOT_ACOUSTIC_CPML,                        &
     &          PML_POTENTIAL_ACOUSTIC_OLD,                             &
     &          PML_POTENTIAL_ACOUSTIC_NEW
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
              INTEGER(KIND=4), INTENT(IN) :: ISPEC_CPML
              REAL(KIND=8), INTENT(IN) :: POTENTIAL_ACOUSTIC(NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: POTENTIAL_DOT_ACOUSTIC(       &
     &NGLOB_AB)
              REAL(KIND=8) :: RMEMORY_POTENTIAL_ACOUSTIC(5,5,5,         &
     &NSPEC_CPML,3)
            END SUBROUTINE PML_COMPUTE_ACCEL_CONTRIBUTION_ACOUSTIC
          END INTERFACE 
        END MODULE PML_COMPUTE_ACCEL_CONTRIBUTION_ACOUSTIC__genmod
