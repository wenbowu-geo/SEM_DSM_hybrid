        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PML_COMPUTE_ACCEL_CONTRIBUTION_ELASTIC__genmod
          INTERFACE 
            SUBROUTINE PML_COMPUTE_ACCEL_CONTRIBUTION_ELASTIC(ISPEC,    &
     &ISPEC_CPML,DISPL,VELOC,RMEMORY_DISPL_ELASTIC)
              USE SPECFEM_PAR, ONLY :                                   &
     &          NGLOB_AB,                                               &
     &          DELTAT,                                                 &
     &          WGLL_CUBE,                                              &
     &          JACOBIAN,                                               &
     &          IBOOL,                                                  &
     &          RHOSTORE,                                               &
     &          CUBED_SPHERE_PROJECTION
              USE PML_PAR, ONLY :                                       &
     &          CPML_REGIONS,                                           &
     &          D_STORE_X,                                              &
     &          D_STORE_Y,                                              &
     &          D_STORE_Z,                                              &
     &          K_STORE_X,                                              &
     &          K_STORE_Y,                                              &
     &          K_STORE_Z,                                              &
     &          ALPHA_STORE_X,                                          &
     &          ALPHA_STORE_Y,                                          &
     &          ALPHA_STORE_Z,                                          &
     &          NSPEC_CPML,                                             &
     &          ACCEL_ELASTIC_CPML,                                     &
     &          PML_DISPL_OLD,                                          &
     &          PML_DISPL_NEW,                                          &
     &          PML_ROTATION_MATRIX
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
              INTEGER(KIND=4), INTENT(IN) :: ISPEC_CPML
              REAL(KIND=8), INTENT(IN) :: DISPL(3,NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: VELOC(3,NGLOB_AB)
              REAL(KIND=8) :: RMEMORY_DISPL_ELASTIC(3,5,5,5,NSPEC_CPML,3&
     &)
            END SUBROUTINE PML_COMPUTE_ACCEL_CONTRIBUTION_ELASTIC
          END INTERFACE 
        END MODULE PML_COMPUTE_ACCEL_CONTRIBUTION_ELASTIC__genmod
