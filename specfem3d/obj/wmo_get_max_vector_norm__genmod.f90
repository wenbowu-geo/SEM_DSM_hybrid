        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:37 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WMO_GET_MAX_VECTOR_NORM__genmod
          INTERFACE 
            SUBROUTINE WMO_GET_MAX_VECTOR_NORM(ISPEC,IGLOB,IA,          &
     &DISPL_ELEMENT,VELOC_ELEMENT,ACCEL_ELEMENT)
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
              INTEGER(KIND=4), INTENT(IN) :: IGLOB
              INTEGER(KIND=4), INTENT(IN) :: IA
              REAL(KIND=8), INTENT(IN) :: DISPL_ELEMENT(3,5,5,5)
              REAL(KIND=8), INTENT(IN) :: VELOC_ELEMENT(3,5,5,5)
              REAL(KIND=8), INTENT(IN) :: ACCEL_ELEMENT(3,5,5,5)
            END SUBROUTINE WMO_GET_MAX_VECTOR_NORM
          END INTERFACE 
        END MODULE WMO_GET_MAX_VECTOR_NORM__genmod
