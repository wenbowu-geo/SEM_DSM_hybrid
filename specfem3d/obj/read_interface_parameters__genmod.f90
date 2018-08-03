        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:49 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_INTERFACE_PARAMETERS__genmod
          INTERFACE 
            SUBROUTINE READ_INTERFACE_PARAMETERS(IUNIT,                 &
     &SUPPRESS_UTM_PROJECTION,INTERFACE_FILE,NPX_INTERFACE,NPY_INTERFACE&
     &,ORIG_X_INTERFACE,ORIG_Y_INTERFACE,SPACING_X_INTERFACE,           &
     &SPACING_Y_INTERFACE,IER)
              INTEGER(KIND=4) :: IUNIT
              LOGICAL(KIND=4) :: SUPPRESS_UTM_PROJECTION
              CHARACTER(LEN=512) :: INTERFACE_FILE
              INTEGER(KIND=4) :: NPX_INTERFACE
              INTEGER(KIND=4) :: NPY_INTERFACE
              REAL(KIND=8) :: ORIG_X_INTERFACE
              REAL(KIND=8) :: ORIG_Y_INTERFACE
              REAL(KIND=8) :: SPACING_X_INTERFACE
              REAL(KIND=8) :: SPACING_Y_INTERFACE
              INTEGER(KIND=4) :: IER
            END SUBROUTINE READ_INTERFACE_PARAMETERS
          END INTERFACE 
        END MODULE READ_INTERFACE_PARAMETERS__genmod
