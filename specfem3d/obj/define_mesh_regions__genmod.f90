        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DEFINE_MESH_REGIONS__genmod
          INTERFACE 
            SUBROUTINE DEFINE_MESH_REGIONS(MYRANK,USE_REGULAR_MESH,     &
     &ISUBREGION,NER,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,IPROC_XI,IPROC_ETA&
     &,NDOUBLINGS,NER_DOUBLINGS,IADDX,IADDY,IADDZ,IX1,IX2,DIX,IY1,IY2,  &
     &DIY,IR1,IR2,DIR,IAX,IAY,IAR)
              INTEGER(KIND=4), INTENT(IN) :: NDOUBLINGS
              INTEGER(KIND=4), INTENT(IN) :: MYRANK
              LOGICAL(KIND=4), INTENT(IN) :: USE_REGULAR_MESH
              INTEGER(KIND=4), INTENT(IN) :: ISUBREGION
              INTEGER(KIND=4), INTENT(IN) :: NER
              INTEGER(KIND=4), INTENT(IN) :: NEX_PER_PROC_XI
              INTEGER(KIND=4), INTENT(IN) :: NEX_PER_PROC_ETA
              INTEGER(KIND=4), INTENT(IN) :: IPROC_XI
              INTEGER(KIND=4), INTENT(IN) :: IPROC_ETA
              INTEGER(KIND=4), INTENT(IN) :: NER_DOUBLINGS(NDOUBLINGS)
              INTEGER(KIND=4), INTENT(INOUT) :: IADDX(8)
              INTEGER(KIND=4), INTENT(INOUT) :: IADDY(8)
              INTEGER(KIND=4), INTENT(INOUT) :: IADDZ(8)
              INTEGER(KIND=4), INTENT(OUT) :: IX1
              INTEGER(KIND=4), INTENT(OUT) :: IX2
              INTEGER(KIND=4), INTENT(OUT) :: DIX
              INTEGER(KIND=4), INTENT(OUT) :: IY1
              INTEGER(KIND=4), INTENT(OUT) :: IY2
              INTEGER(KIND=4), INTENT(OUT) :: DIY
              INTEGER(KIND=4), INTENT(OUT) :: IR1
              INTEGER(KIND=4), INTENT(OUT) :: IR2
              INTEGER(KIND=4), INTENT(OUT) :: DIR
              INTEGER(KIND=4), INTENT(OUT) :: IAX
              INTEGER(KIND=4), INTENT(OUT) :: IAY
              INTEGER(KIND=4), INTENT(OUT) :: IAR
            END SUBROUTINE DEFINE_MESH_REGIONS
          END INTERFACE 
        END MODULE DEFINE_MESH_REGIONS__genmod
