        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:43 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_PARAMETERS__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_PARAMETERS(NER,NEX_XI,NEX_ETA,NPROC_XI,  &
     &NPROC_ETA,NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NSPEC_AB,        &
     &NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_A_ETA,NSPEC2D_B_ETA,            &
     &NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,         &
     &NSPEC2D_TOP,NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB,   &
     &USE_REGULAR_MESH,NDOUBLINGS,NER_DOUBLINGS)
              INTEGER(KIND=4), INTENT(IN) :: NDOUBLINGS
              INTEGER(KIND=4), INTENT(IN) :: NER
              INTEGER(KIND=4), INTENT(IN) :: NEX_XI
              INTEGER(KIND=4), INTENT(IN) :: NEX_ETA
              INTEGER(KIND=4), INTENT(IN) :: NPROC_XI
              INTEGER(KIND=4), INTENT(IN) :: NPROC_ETA
              INTEGER(KIND=4), INTENT(OUT) :: NPROC
              INTEGER(KIND=4), INTENT(OUT) :: NEX_PER_PROC_XI
              INTEGER(KIND=4), INTENT(OUT) :: NEX_PER_PROC_ETA
              INTEGER(KIND=4), INTENT(OUT) :: NSPEC_AB
              INTEGER(KIND=4), INTENT(OUT) :: NSPEC2D_A_XI
              INTEGER(KIND=4), INTENT(OUT) :: NSPEC2D_B_XI
              INTEGER(KIND=4), INTENT(OUT) :: NSPEC2D_A_ETA
              INTEGER(KIND=4), INTENT(OUT) :: NSPEC2D_B_ETA
              INTEGER(KIND=4), INTENT(OUT) :: NSPEC2DMAX_XMIN_XMAX
              INTEGER(KIND=4), INTENT(OUT) :: NSPEC2DMAX_YMIN_YMAX
              INTEGER(KIND=4), INTENT(OUT) :: NSPEC2D_BOTTOM
              INTEGER(KIND=4), INTENT(OUT) :: NSPEC2D_TOP
              INTEGER(KIND=4), INTENT(OUT) :: NPOIN2DMAX_XMIN_XMAX
              INTEGER(KIND=4), INTENT(OUT) :: NPOIN2DMAX_YMIN_YMAX
              INTEGER(KIND=4), INTENT(OUT) :: NGLOB_AB
              LOGICAL(KIND=4), INTENT(IN) :: USE_REGULAR_MESH
              INTEGER(KIND=4), INTENT(IN) :: NER_DOUBLINGS(NDOUBLINGS)
            END SUBROUTINE COMPUTE_PARAMETERS
          END INTERFACE 
        END MODULE COMPUTE_PARAMETERS__genmod