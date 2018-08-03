        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:30 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DS_SET_CROSS_SECTION_FLAGS__genmod
          INTERFACE 
            SUBROUTINE DS_SET_CROSS_SECTION_FLAGS(NSPEC,                &
     &ISPEC_IS_SURFACE_EXTERNAL_MESH,NGLOB,                             &
     &IGLOB_IS_SURFACE_EXTERNAL_MESH,I,J,K,ISPEC,IBOOL,                 &
     &VALENCE_EXTERNAL_MESH,COUNTVAL)
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: NSPEC
              LOGICAL(KIND=4) :: ISPEC_IS_SURFACE_EXTERNAL_MESH(NSPEC)
              LOGICAL(KIND=4) :: IGLOB_IS_SURFACE_EXTERNAL_MESH(NGLOB)
              INTEGER(KIND=4) :: I
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: ISPEC
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              INTEGER(KIND=4) :: VALENCE_EXTERNAL_MESH(NGLOB)
              INTEGER(KIND=4) :: COUNTVAL
            END SUBROUTINE DS_SET_CROSS_SECTION_FLAGS
          END INTERFACE 
        END MODULE DS_SET_CROSS_SECTION_FLAGS__genmod
