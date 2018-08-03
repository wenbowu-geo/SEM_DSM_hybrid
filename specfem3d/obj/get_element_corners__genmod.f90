        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:33 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ELEMENT_CORNERS__genmod
          INTERFACE 
            SUBROUTINE GET_ELEMENT_CORNERS(ISPEC,IFACE_REF,XCOORD,YCOORD&
     &,ZCOORD,IGLOB_CORNERS_REF,IBOOL,NSPEC,NGLOB,XSTORE_DUMMY,         &
     &YSTORE_DUMMY,ZSTORE_DUMMY,IFACE_ALL_CORNER_IJK)
              INTEGER(KIND=4), INTENT(IN) :: NGLOB
              INTEGER(KIND=4), INTENT(IN) :: NSPEC
              INTEGER(KIND=4), INTENT(IN) :: ISPEC
              INTEGER(KIND=4), INTENT(IN) :: IFACE_REF
              REAL(KIND=8), INTENT(OUT) :: XCOORD(4)
              REAL(KIND=8), INTENT(OUT) :: YCOORD(4)
              REAL(KIND=8), INTENT(OUT) :: ZCOORD(4)
              INTEGER(KIND=4), INTENT(OUT) :: IGLOB_CORNERS_REF(4)
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              REAL(KIND=8) :: XSTORE_DUMMY(NGLOB)
              REAL(KIND=8) :: YSTORE_DUMMY(NGLOB)
              REAL(KIND=8) :: ZSTORE_DUMMY(NGLOB)
              INTEGER(KIND=4) :: IFACE_ALL_CORNER_IJK(3,4,6)
            END SUBROUTINE GET_ELEMENT_CORNERS
          END INTERFACE 
        END MODULE GET_ELEMENT_CORNERS__genmod
