        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:20 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PML_SET_LOCAL_DAMPINGCOEFF__genmod
          INTERFACE 
            SUBROUTINE PML_SET_LOCAL_DAMPINGCOEFF(MYRANK,XSTORE,YSTORE, &
     &ZSTORE)
              USE GENERATE_DATABASES_PAR, ONLY :                        &
     &          IBOOL,                                                  &
     &          NGLOB_AB,                                               &
     &          D_STORE_X,                                              &
     &          D_STORE_Y,                                              &
     &          D_STORE_Z,                                              &
     &          K_STORE_X,                                              &
     &          K_STORE_Y,                                              &
     &          K_STORE_Z,                                              &
     &          ALPHA_STORE_X,                                          &
     &          ALPHA_STORE_Y,                                          &
     &          ALPHA_STORE_Z,                                          &
     &          CPML_TO_SPEC,                                           &
     &          CPML_WIDTH_X,                                           &
     &          CPML_WIDTH_Y,                                           &
     &          CPML_WIDTH_Z,                                           &
     &          MIN_DISTANCE_BETWEEN_CPML_PARAMETER,                    &
     &          NPOWER,                                                 &
     &          CUSTOM_REAL,                                            &
     &          SIZE_REAL,                                              &
     &          NGLLX,                                                  &
     &          NGLLY,                                                  &
     &          NGLLZ,                                                  &
     &          NSPEC_CPML,                                             &
     &          PML_INSTEAD_OF_FREE_SURFACE,                            &
     &          IMAIN,                                                  &
     &          CPML_REGIONS,                                           &
     &          F0_FOR_PML,                                             &
     &          PI,                                                     &
     &          CPML_X_ONLY,                                            &
     &          CPML_Y_ONLY,                                            &
     &          CPML_Z_ONLY,                                            &
     &          CPML_XY_ONLY,                                           &
     &          CPML_XZ_ONLY,                                           &
     &          CPML_YZ_ONLY,                                           &
     &          CPML_XYZ,                                               &
     &          SIMULATION_TYPE,                                        &
     &          SAVE_FORWARD,                                           &
     &          IS_CPML,                                                &
     &          MASK_IBOOL_INTERIOR_DOMAIN,                             &
     &          NGLOB_INTERFACE_PML_ACOUSTIC,                           &
     &          POINTS_INTERFACE_PML_ACOUSTIC,                          &
     &          NGLOB_INTERFACE_PML_ELASTIC,                            &
     &          POINTS_INTERFACE_PML_ELASTIC,                           &
     &          ZERO,                                                   &
     &          ONE,                                                    &
     &          TWO,                                                    &
     &          HUGEVAL,                                                &
     &          COUPLING_TYPE, ONLY :                                   &
     &          NSPEC =>NSPEC_AB
              INTEGER(KIND=4), INTENT(IN) :: MYRANK
              REAL(KIND=8), INTENT(IN) :: XSTORE(NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: YSTORE(NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: ZSTORE(NGLOB_AB)
            END SUBROUTINE PML_SET_LOCAL_DAMPINGCOEFF
          END INTERFACE 
        END MODULE PML_SET_LOCAL_DAMPINGCOEFF__genmod
