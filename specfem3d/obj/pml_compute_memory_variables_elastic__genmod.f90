        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:19 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PML_COMPUTE_MEMORY_VARIABLES_ELASTIC__genmod
          INTERFACE 
            SUBROUTINE PML_COMPUTE_MEMORY_VARIABLES_ELASTIC(ISPEC,      &
     &ISPEC_CPML,TEMPX1,TEMPY1,TEMPZ1,TEMPX2,TEMPY2,TEMPZ2,TEMPX3,TEMPY3&
     &,TEMPZ3,RMEMORY_DUX_DXL_X,RMEMORY_DUY_DYL_X,RMEMORY_DUZ_DZL_X,    &
     &RMEMORY_DUX_DYL_X,RMEMORY_DUX_DZL_X,RMEMORY_DUZ_DXL_X,            &
     &RMEMORY_DUY_DXL_X,RMEMORY_DUX_DXL_Y,RMEMORY_DUZ_DZL_Y,            &
     &RMEMORY_DUY_DYL_Y,RMEMORY_DUY_DXL_Y,RMEMORY_DUY_DZL_Y,            &
     &RMEMORY_DUZ_DYL_Y,RMEMORY_DUX_DYL_Y,RMEMORY_DUX_DXL_Z,            &
     &RMEMORY_DUY_DYL_Z,RMEMORY_DUZ_DZL_Z,RMEMORY_DUZ_DXL_Z,            &
     &RMEMORY_DUZ_DYL_Z,RMEMORY_DUY_DZL_Z,RMEMORY_DUX_DZL_Z)
              USE PML_PAR, ONLY :                                       &
     &          NSPEC_CPML,                                             &
     &          CPML_REGIONS,                                           &
     &          K_STORE_X,                                              &
     &          K_STORE_Y,                                              &
     &          K_STORE_Z,                                              &
     &          D_STORE_X,                                              &
     &          D_STORE_Y,                                              &
     &          D_STORE_Z,                                              &
     &          ALPHA_STORE_X,                                          &
     &          ALPHA_STORE_Y,                                          &
     &          ALPHA_STORE_Z,                                          &
     &          PML_DUX_DXL,                                            &
     &          PML_DUX_DYL,                                            &
     &          PML_DUX_DZL,                                            &
     &          PML_DUY_DXL,                                            &
     &          PML_DUY_DYL,                                            &
     &          PML_DUY_DZL,                                            &
     &          PML_DUZ_DXL,                                            &
     &          PML_DUZ_DYL,                                            &
     &          PML_DUZ_DZL,                                            &
     &          PML_DUX_DXL_OLD,                                        &
     &          PML_DUX_DYL_OLD,                                        &
     &          PML_DUX_DZL_OLD,                                        &
     &          PML_DUY_DXL_OLD,                                        &
     &          PML_DUY_DYL_OLD,                                        &
     &          PML_DUY_DZL_OLD,                                        &
     &          PML_DUZ_DXL_OLD,                                        &
     &          PML_DUZ_DYL_OLD,                                        &
     &          PML_DUZ_DZL_OLD,                                        &
     &          PML_DUX_DXL_NEW,                                        &
     &          PML_DUX_DYL_NEW,                                        &
     &          PML_DUX_DZL_NEW,                                        &
     &          PML_DUY_DXL_NEW,                                        &
     &          PML_DUY_DYL_NEW,                                        &
     &          PML_DUY_DZL_NEW,                                        &
     &          PML_DUZ_DXL_NEW,                                        &
     &          PML_DUZ_DYL_NEW,                                        &
     &          PML_DUZ_DZL_NEW,                                        &
     &          PML_ROTATION_MATRIX
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
              INTEGER(KIND=4), INTENT(IN) :: ISPEC_CPML
              REAL(KIND=8), INTENT(OUT) :: TEMPX1(5,5,5)
              REAL(KIND=8), INTENT(OUT) :: TEMPY1(5,5,5)
              REAL(KIND=8), INTENT(OUT) :: TEMPZ1(5,5,5)
              REAL(KIND=8), INTENT(OUT) :: TEMPX2(5,5,5)
              REAL(KIND=8), INTENT(OUT) :: TEMPY2(5,5,5)
              REAL(KIND=8), INTENT(OUT) :: TEMPZ2(5,5,5)
              REAL(KIND=8), INTENT(OUT) :: TEMPX3(5,5,5)
              REAL(KIND=8), INTENT(OUT) :: TEMPY3(5,5,5)
              REAL(KIND=8), INTENT(OUT) :: TEMPZ3(5,5,5)
              REAL(KIND=8) :: RMEMORY_DUX_DXL_X(5,5,5,NSPEC_CPML,3)
              REAL(KIND=8) :: RMEMORY_DUY_DYL_X(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUZ_DZL_X(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUX_DYL_X(5,5,5,NSPEC_CPML,3)
              REAL(KIND=8) :: RMEMORY_DUX_DZL_X(5,5,5,NSPEC_CPML,3)
              REAL(KIND=8) :: RMEMORY_DUZ_DXL_X(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUY_DXL_X(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUX_DXL_Y(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUZ_DZL_Y(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUY_DYL_Y(5,5,5,NSPEC_CPML,3)
              REAL(KIND=8) :: RMEMORY_DUY_DXL_Y(5,5,5,NSPEC_CPML,3)
              REAL(KIND=8) :: RMEMORY_DUY_DZL_Y(5,5,5,NSPEC_CPML,3)
              REAL(KIND=8) :: RMEMORY_DUZ_DYL_Y(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUX_DYL_Y(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUX_DXL_Z(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUY_DYL_Z(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUZ_DZL_Z(5,5,5,NSPEC_CPML,3)
              REAL(KIND=8) :: RMEMORY_DUZ_DXL_Z(5,5,5,NSPEC_CPML,3)
              REAL(KIND=8) :: RMEMORY_DUZ_DYL_Z(5,5,5,NSPEC_CPML,3)
              REAL(KIND=8) :: RMEMORY_DUY_DZL_Z(5,5,5,NSPEC_CPML)
              REAL(KIND=8) :: RMEMORY_DUX_DZL_Z(5,5,5,NSPEC_CPML)
            END SUBROUTINE PML_COMPUTE_MEMORY_VARIABLES_ELASTIC
          END INTERFACE 
        END MODULE PML_COMPUTE_MEMORY_VARIABLES_ELASTIC__genmod
