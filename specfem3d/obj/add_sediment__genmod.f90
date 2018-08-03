        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:56 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ADD_SEDIMENT__genmod
          INTERFACE 
            SUBROUTINE ADD_SEDIMENT(MYRANK,X_EVAL,Y_EVAL,Z_EVAL,        &
     &FLAG_MEDIA,FLUID_THISPOINT,VP_FINAL,VS_FINAL,RHO_FINAL)
              INTEGER(KIND=4), INTENT(IN) :: MYRANK
              REAL(KIND=8), INTENT(IN) :: X_EVAL
              REAL(KIND=8), INTENT(IN) :: Y_EVAL
              REAL(KIND=8), INTENT(IN) :: Z_EVAL
              INTEGER(KIND=4), INTENT(IN) :: FLAG_MEDIA
              INTEGER(KIND=4), INTENT(IN) :: FLUID_THISPOINT
              REAL(KIND=8), INTENT(OUT) :: VP_FINAL
              REAL(KIND=8), INTENT(OUT) :: VS_FINAL
              REAL(KIND=8), INTENT(OUT) :: RHO_FINAL
            END SUBROUTINE ADD_SEDIMENT
          END INTERFACE 
        END MODULE ADD_SEDIMENT__genmod
