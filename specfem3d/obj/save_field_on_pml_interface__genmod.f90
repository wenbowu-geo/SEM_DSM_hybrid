        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SAVE_FIELD_ON_PML_INTERFACE__genmod
          INTERFACE 
            SUBROUTINE SAVE_FIELD_ON_PML_INTERFACE(DISPL,VELOC,ACCEL,   &
     &NGLOB_INTERFACE_PML_ELASTIC,B_PML_FIELD,B_RECLEN_PML_FIELD)
              USE SPECFEM_PAR, ONLY :                                   &
     &          NGLOB_AB,                                               &
     &          IT
              INTEGER(KIND=4), INTENT(IN) :: NGLOB_INTERFACE_PML_ELASTIC
              REAL(KIND=8), INTENT(IN) :: DISPL(3,NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: VELOC(3,NGLOB_AB)
              REAL(KIND=8), INTENT(IN) :: ACCEL(3,NGLOB_AB)
              REAL(KIND=8) :: B_PML_FIELD(9,NGLOB_INTERFACE_PML_ELASTIC)
              INTEGER(KIND=4), INTENT(IN) :: B_RECLEN_PML_FIELD
            END SUBROUTINE SAVE_FIELD_ON_PML_INTERFACE
          END INTERFACE 
        END MODULE SAVE_FIELD_ON_PML_INTERFACE__genmod
