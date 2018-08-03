        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:19 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PML_COMPUTE_MEMORY_VARIABLES_ELASTIC_ACOUSTIC__genmod
          INTERFACE 
            SUBROUTINE PML_COMPUTE_MEMORY_VARIABLES_ELASTIC_ACOUSTIC(   &
     &ISPEC_CPML,IFACE,IGLOB,I,J,K,PRESSURE,POTENTIAL_ACOUSTIC,         &
     &POTENTIAL_DOT_ACOUSTIC,POTENTIAL_DOT_DOT_ACOUSTIC,                &
     &NUM_COUPLING_AC_EL_FACES,RMEMORY_COUPLING_EL_AC_POTENTIAL)
              USE SPECFEM_PAR, ONLY :                                   &
     &          NGLOB_AB,                                               &
     &          DELTAT
              INTEGER(KIND=4), INTENT(IN) :: NUM_COUPLING_AC_EL_FACES
              INTEGER(KIND=4), INTENT(IN) :: ISPEC_CPML
              INTEGER(KIND=4), INTENT(IN) :: IFACE
              INTEGER(KIND=4), INTENT(IN) :: IGLOB
              INTEGER(KIND=4) :: I
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: K
              REAL(KIND=8), INTENT(OUT) :: PRESSURE
              REAL(KIND=8), INTENT(IN) :: POTENTIAL_ACOUSTIC(NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: POTENTIAL_DOT_ACOUSTIC(       &
     &NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: POTENTIAL_DOT_DOT_ACOUSTIC(   &
     &NGLOB_AB)
              REAL(KIND=8) :: RMEMORY_COUPLING_EL_AC_POTENTIAL(5,5,5,   &
     &NUM_COUPLING_AC_EL_FACES,2)
            END SUBROUTINE PML_COMPUTE_MEMORY_VARIABLES_ELASTIC_ACOUSTIC
          END INTERFACE 
        END MODULE PML_COMPUTE_MEMORY_VARIABLES_ELASTIC_ACOUSTIC__genmod
