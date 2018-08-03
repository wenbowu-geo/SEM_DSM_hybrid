        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:19 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PML_COMPUTE_MEMORY_VARIABLES_ACOUSTIC__genmod
          INTERFACE 
            SUBROUTINE PML_COMPUTE_MEMORY_VARIABLES_ACOUSTIC(ISPEC,     &
     &ISPEC_CPML,TEMP1,TEMP2,TEMP3,RMEMORY_DPOTENTIAL_DXL,              &
     &RMEMORY_DPOTENTIAL_DYL,RMEMORY_DPOTENTIAL_DZL,PML_DPOTENTIAL_DXL, &
     &PML_DPOTENTIAL_DYL,PML_DPOTENTIAL_DZL,PML_DPOTENTIAL_DXL_OLD,     &
     &PML_DPOTENTIAL_DYL_OLD,PML_DPOTENTIAL_DZL_OLD,                    &
     &PML_DPOTENTIAL_DXL_NEW,PML_DPOTENTIAL_DYL_NEW,                    &
     &PML_DPOTENTIAL_DZL_NEW)
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
     &          ALPHA_STORE_Z
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
              INTEGER(KIND=4), INTENT(IN) :: ISPEC_CPML
              REAL(KIND=8), INTENT(OUT) :: TEMP1(5,5,5)
              REAL(KIND=8), INTENT(OUT) :: TEMP2(5,5,5)
              REAL(KIND=8), INTENT(OUT) :: TEMP3(5,5,5)
              REAL(KIND=8), INTENT(INOUT) :: RMEMORY_DPOTENTIAL_DXL(5,5,&
     &5,NSPEC_CPML,3)
              REAL(KIND=8), INTENT(INOUT) :: RMEMORY_DPOTENTIAL_DYL(5,5,&
     &5,NSPEC_CPML,3)
              REAL(KIND=8), INTENT(INOUT) :: RMEMORY_DPOTENTIAL_DZL(5,5,&
     &5,NSPEC_CPML,3)
              REAL(KIND=8), INTENT(IN) :: PML_DPOTENTIAL_DXL(5,5,5)
              REAL(KIND=8), INTENT(IN) :: PML_DPOTENTIAL_DYL(5,5,5)
              REAL(KIND=8), INTENT(IN) :: PML_DPOTENTIAL_DZL(5,5,5)
              REAL(KIND=8), INTENT(IN) :: PML_DPOTENTIAL_DXL_OLD(5,5,5)
              REAL(KIND=8), INTENT(IN) :: PML_DPOTENTIAL_DYL_OLD(5,5,5)
              REAL(KIND=8), INTENT(IN) :: PML_DPOTENTIAL_DZL_OLD(5,5,5)
              REAL(KIND=8), INTENT(IN) :: PML_DPOTENTIAL_DXL_NEW(5,5,5)
              REAL(KIND=8), INTENT(IN) :: PML_DPOTENTIAL_DYL_NEW(5,5,5)
              REAL(KIND=8), INTENT(IN) :: PML_DPOTENTIAL_DZL_NEW(5,5,5)
            END SUBROUTINE PML_COMPUTE_MEMORY_VARIABLES_ACOUSTIC
          END INTERFACE 
        END MODULE PML_COMPUTE_MEMORY_VARIABLES_ACOUSTIC__genmod
