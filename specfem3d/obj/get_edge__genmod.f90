        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:34 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_EDGE__genmod
          INTERFACE 
            SUBROUTINE GET_EDGE(N,ITYPE,E1,E2,E3,E4,IXMIN,IXMAX,IYMIN,  &
     &IYMAX,IZMIN,IZMAX)
              INTEGER(KIND=4), INTENT(IN) :: N(8)
              INTEGER(KIND=4), INTENT(IN) :: ITYPE
              INTEGER(KIND=4), INTENT(IN) :: E1
              INTEGER(KIND=4), INTENT(IN) :: E2
              INTEGER(KIND=4), INTENT(IN) :: E3
              INTEGER(KIND=4), INTENT(IN) :: E4
              INTEGER(KIND=4), INTENT(OUT) :: IXMIN
              INTEGER(KIND=4), INTENT(OUT) :: IXMAX
              INTEGER(KIND=4), INTENT(OUT) :: IYMIN
              INTEGER(KIND=4), INTENT(OUT) :: IYMAX
              INTEGER(KIND=4), INTENT(OUT) :: IZMIN
              INTEGER(KIND=4), INTENT(OUT) :: IZMAX
            END SUBROUTINE GET_EDGE
          END INTERFACE 
        END MODULE GET_EDGE__genmod
