        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DEFINE_MODEL_REGIONS__genmod
          INTERFACE 
            SUBROUTINE DEFINE_MODEL_REGIONS(NEX_PER_PROC_XI,            &
     &NEX_PER_PROC_ETA,IPROC_XI,IPROC_ETA,ISUBREGION,NBSUBREGIONS,      &
     &SUBREGIONS,IADDX,IADDY,IADDZ,IX1,IX2,DIX,IY1,IY2,DIY,IR1,IR2,DIR, &
     &IAX,IAY,IAR,NUM_MATERIAL)
              INTEGER(KIND=4), INTENT(IN) :: NBSUBREGIONS
              INTEGER(KIND=4), INTENT(IN) :: NEX_PER_PROC_XI
              INTEGER(KIND=4), INTENT(IN) :: NEX_PER_PROC_ETA
              INTEGER(KIND=4), INTENT(IN) :: IPROC_XI
              INTEGER(KIND=4), INTENT(IN) :: IPROC_ETA
              INTEGER(KIND=4), INTENT(IN) :: ISUBREGION
              INTEGER(KIND=4), INTENT(IN) :: SUBREGIONS(NBSUBREGIONS,7)
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
              INTEGER(KIND=4), INTENT(OUT) :: NUM_MATERIAL
            END SUBROUTINE DEFINE_MODEL_REGIONS
          END INTERFACE 
        END MODULE DEFINE_MODEL_REGIONS__genmod
