        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:38 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SORT_ARRAY_COORDINATES__genmod
          INTERFACE 
            SUBROUTINE SORT_ARRAY_COORDINATES(NPOINTOT,X,Y,Z,IBOOL,IGLOB&
     &,LOCVAL,IFSEG,NGLOB,NINSEG,XTOL)
              INTEGER(KIND=4), INTENT(IN) :: NPOINTOT
              REAL(KIND=8), INTENT(INOUT) :: X(NPOINTOT)
              REAL(KIND=8), INTENT(INOUT) :: Y(NPOINTOT)
              REAL(KIND=8), INTENT(INOUT) :: Z(NPOINTOT)
              INTEGER(KIND=4), INTENT(INOUT) :: IBOOL(NPOINTOT)
              INTEGER(KIND=4), INTENT(OUT) :: IGLOB(NPOINTOT)
              INTEGER(KIND=4), INTENT(OUT) :: LOCVAL(NPOINTOT)
              LOGICAL(KIND=4), INTENT(OUT) :: IFSEG(NPOINTOT)
              INTEGER(KIND=4), INTENT(OUT) :: NGLOB
              INTEGER(KIND=4), INTENT(OUT) :: NINSEG(NPOINTOT)
              REAL(KIND=8), INTENT(IN) :: XTOL
            END SUBROUTINE SORT_ARRAY_COORDINATES
          END INTERFACE 
        END MODULE SORT_ARRAY_COORDINATES__genmod
