        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:59 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CRM_EXT_ALLOCATE_ARRAYS__genmod
          INTERFACE 
            SUBROUTINE CRM_EXT_ALLOCATE_ARRAYS(NSPEC,LOCAL_PATH,MYRANK, &
     &NSPEC2D_XMIN,NSPEC2D_XMAX,NSPEC2D_YMIN,NSPEC2D_YMAX,NSPEC2D_BOTTOM&
     &,NSPEC2D_TOP,ANISOTROPY)
              INTEGER(KIND=4) :: NSPEC
              CHARACTER(LEN=512) :: LOCAL_PATH
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: NSPEC2D_XMIN
              INTEGER(KIND=4) :: NSPEC2D_XMAX
              INTEGER(KIND=4) :: NSPEC2D_YMIN
              INTEGER(KIND=4) :: NSPEC2D_YMAX
              INTEGER(KIND=4) :: NSPEC2D_BOTTOM
              INTEGER(KIND=4) :: NSPEC2D_TOP
              LOGICAL(KIND=4) :: ANISOTROPY
            END SUBROUTINE CRM_EXT_ALLOCATE_ARRAYS
          END INTERFACE 
        END MODULE CRM_EXT_ALLOCATE_ARRAYS__genmod
