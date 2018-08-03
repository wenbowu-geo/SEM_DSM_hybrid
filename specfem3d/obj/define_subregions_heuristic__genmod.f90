        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DEFINE_SUBREGIONS_HEURISTIC__genmod
          INTERFACE 
            SUBROUTINE DEFINE_SUBREGIONS_HEURISTIC(MYRANK,ISUBREGION,   &
     &IADDX,IADDY,IADDZ,IX1,IX2,DIX,IY1,IY2,DIY,IR1,IR2,DIR,IAX,IAY,IAR,&
     &ITYPE_ELEMENT,NPX,NPY,NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,&
     &NER_BASEMENT_SEDIM)
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: ISUBREGION
              INTEGER(KIND=4) :: IADDX(8)
              INTEGER(KIND=4) :: IADDY(8)
              INTEGER(KIND=4) :: IADDZ(8)
              INTEGER(KIND=4) :: IX1
              INTEGER(KIND=4) :: IX2
              INTEGER(KIND=4) :: DIX
              INTEGER(KIND=4) :: IY1
              INTEGER(KIND=4) :: IY2
              INTEGER(KIND=4) :: DIY
              INTEGER(KIND=4) :: IR1
              INTEGER(KIND=4) :: IR2
              INTEGER(KIND=4) :: DIR
              INTEGER(KIND=4) :: IAX
              INTEGER(KIND=4) :: IAY
              INTEGER(KIND=4) :: IAR
              INTEGER(KIND=4) :: ITYPE_ELEMENT
              INTEGER(KIND=4) :: NPX
              INTEGER(KIND=4) :: NPY
              INTEGER(KIND=4) :: NER_BOTTOM_MOHO
              INTEGER(KIND=4) :: NER_MOHO_16
              INTEGER(KIND=4) :: NER_16_BASEMENT
              INTEGER(KIND=4) :: NER_BASEMENT_SEDIM
            END SUBROUTINE DEFINE_SUBREGIONS_HEURISTIC
          END INTERFACE 
        END MODULE DEFINE_SUBREGIONS_HEURISTIC__genmod
