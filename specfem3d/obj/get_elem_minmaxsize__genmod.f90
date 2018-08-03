        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:28 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ELEM_MINMAXSIZE__genmod
          INTERFACE 
            SUBROUTINE GET_ELEM_MINMAXSIZE(ELEMSIZE_MIN,ELEMSIZE_MAX,   &
     &ISPEC,NSPEC_AB,NGLOB_AB,IBOOL,XSTORE,YSTORE,ZSTORE)
              INTEGER(KIND=4) :: NGLOB_AB
              INTEGER(KIND=4) :: NSPEC_AB
              REAL(KIND=8) :: ELEMSIZE_MIN
              REAL(KIND=8) :: ELEMSIZE_MAX
              INTEGER(KIND=4) :: ISPEC
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: XSTORE(NGLOB_AB)
              REAL(KIND=8) :: YSTORE(NGLOB_AB)
              REAL(KIND=8) :: ZSTORE(NGLOB_AB)
            END SUBROUTINE GET_ELEM_MINMAXSIZE
          END INTERFACE 
        END MODULE GET_ELEM_MINMAXSIZE__genmod
