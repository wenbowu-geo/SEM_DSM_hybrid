        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:33 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ELEMENT_FACE_GLL_INDICES__genmod
          INTERFACE 
            SUBROUTINE GET_ELEMENT_FACE_GLL_INDICES(IFACE,IJK_FACE,NGLLA&
     &,NGLLB)
              INTEGER(KIND=4) :: NGLLB
              INTEGER(KIND=4) :: NGLLA
              INTEGER(KIND=4) :: IFACE
              INTEGER(KIND=4) :: IJK_FACE(3,NGLLA,NGLLB)
            END SUBROUTINE GET_ELEMENT_FACE_GLL_INDICES
          END INTERFACE 
        END MODULE GET_ELEMENT_FACE_GLL_INDICES__genmod
