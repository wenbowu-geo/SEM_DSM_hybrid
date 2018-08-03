        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:19 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PML_COMPUTE_MEMORY_VARIABLES_ACOUSTIC_ELASTIC__genmod
          INTERFACE 
            SUBROUTINE PML_COMPUTE_MEMORY_VARIABLES_ACOUSTIC_ELASTIC(   &
     &ISPEC_CPML,IFACE,IGLOB,I,J,K,DISPL_X,DISPL_Y,DISPL_Z,DISPL,       &
     &NUM_COUPLING_AC_EL_FACES,RMEMORY_COUPLING_AC_EL_DISPL)
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
              REAL(KIND=8), INTENT(OUT) :: DISPL_X
              REAL(KIND=8), INTENT(OUT) :: DISPL_Y
              REAL(KIND=8), INTENT(OUT) :: DISPL_Z
              REAL(KIND=8), INTENT(IN) :: DISPL(3,NGLOB_AB)
              REAL(KIND=8), INTENT(INOUT) ::                            &
     &RMEMORY_COUPLING_AC_EL_DISPL(3,5,5,5,NUM_COUPLING_AC_EL_FACES,2)
            END SUBROUTINE PML_COMPUTE_MEMORY_VARIABLES_ACOUSTIC_ELASTIC
          END INTERFACE 
        END MODULE PML_COMPUTE_MEMORY_VARIABLES_ACOUSTIC_ELASTIC__genmod
