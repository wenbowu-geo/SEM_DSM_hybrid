        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:30 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DS_SET_PLANE_FLAGS__genmod
          INTERFACE 
            SUBROUTINE DS_SET_PLANE_FLAGS(IFACE,ISPEC,NSPEC,            &
     &ISPEC_IS_SURFACE_EXTERNAL_MESH,NGLOB,                             &
     &IGLOB_IS_SURFACE_EXTERNAL_MESH,IBOOL,VALENCE_EXTERNAL_MESH)
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: IFACE
              INTEGER(KIND=4) :: ISPEC
              LOGICAL(KIND=4) :: ISPEC_IS_SURFACE_EXTERNAL_MESH(NSPEC)
              LOGICAL(KIND=4) :: IGLOB_IS_SURFACE_EXTERNAL_MESH(NGLOB)
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              INTEGER(KIND=4) :: VALENCE_EXTERNAL_MESH(NGLOB)
            END SUBROUTINE DS_SET_PLANE_FLAGS
          END INTERFACE 
        END MODULE DS_SET_PLANE_FLAGS__genmod
