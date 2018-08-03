        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:37 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WMO_GET_VEL_VECTOR__genmod
          INTERFACE 
            SUBROUTINE WMO_GET_VEL_VECTOR(ISPEC,IGLOB,IA,VAL_ELEMENT)
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
              INTEGER(KIND=4), INTENT(IN) :: IGLOB
              INTEGER(KIND=4), INTENT(IN) :: IA
              REAL(KIND=8), INTENT(IN) :: VAL_ELEMENT(3,5,5,5)
            END SUBROUTINE WMO_GET_VEL_VECTOR
          END INTERFACE 
        END MODULE WMO_GET_VEL_VECTOR__genmod
