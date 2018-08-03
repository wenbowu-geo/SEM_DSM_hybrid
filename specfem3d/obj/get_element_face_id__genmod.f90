        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:33 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ELEMENT_FACE_ID__genmod
          INTERFACE 
            SUBROUTINE GET_ELEMENT_FACE_ID(ISPEC,XCOORD,YCOORD,ZCOORD,  &
     &IBOOL,NSPEC,NGLOB,XSTORE_DUMMY,YSTORE_DUMMY,ZSTORE_DUMMY,IFACE_ID)
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: ISPEC
              REAL(KIND=8) :: XCOORD(4)
              REAL(KIND=8) :: YCOORD(4)
              REAL(KIND=8) :: ZCOORD(4)
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              REAL(KIND=8) :: XSTORE_DUMMY(NGLOB)
              REAL(KIND=8) :: YSTORE_DUMMY(NGLOB)
              REAL(KIND=8) :: ZSTORE_DUMMY(NGLOB)
              INTEGER(KIND=4) :: IFACE_ID
            END SUBROUTINE GET_ELEMENT_FACE_ID
          END INTERFACE 
        END MODULE GET_ELEMENT_FACE_ID__genmod
