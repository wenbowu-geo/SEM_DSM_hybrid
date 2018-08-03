        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:51 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_GLOBAL__genmod
          INTERFACE 
            SUBROUTINE GET_GLOBAL(NPOINTOT,XP,YP,ZP,IGLOB,LOCVAL,IFSEG, &
     &NGLOB,UTM_X_MIN,UTM_X_MAX)
              INTEGER(KIND=4) :: NPOINTOT
              REAL(KIND=8) :: XP(NPOINTOT)
              REAL(KIND=8) :: YP(NPOINTOT)
              REAL(KIND=8) :: ZP(NPOINTOT)
              INTEGER(KIND=4) :: IGLOB(NPOINTOT)
              INTEGER(KIND=4) :: LOCVAL(NPOINTOT)
              LOGICAL(KIND=4) :: IFSEG(NPOINTOT)
              INTEGER(KIND=4) :: NGLOB
              REAL(KIND=8) :: UTM_X_MIN
              REAL(KIND=8) :: UTM_X_MAX
            END SUBROUTINE GET_GLOBAL
          END INTERFACE 
        END MODULE GET_GLOBAL__genmod
