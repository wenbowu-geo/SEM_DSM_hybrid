        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:01 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ABSORBING_BOUNDARY__genmod
          INTERFACE 
            SUBROUTINE GET_ABSORBING_BOUNDARY(MYRANK,NSPEC,IBOOL,       &
     &NODES_COORDS_EXT_MESH,NNODES_EXT_MESH,IBELM_XMIN,IBELM_XMAX,      &
     &IBELM_YMIN,IBELM_YMAX,IBELM_BOTTOM,IBELM_TOP,NODES_IBELM_XMIN,    &
     &NODES_IBELM_XMAX,NODES_IBELM_YMIN,NODES_IBELM_YMAX,               &
     &NODES_IBELM_BOTTOM,NODES_IBELM_TOP,NSPEC2D_XMIN,NSPEC2D_XMAX,     &
     &NSPEC2D_YMIN,NSPEC2D_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP)
              USE GENERATE_DATABASES_PAR, ONLY :                        &
     &          STACEY_INSTEAD_OF_FREE_SURFACE,                         &
     &          PML_INSTEAD_OF_FREE_SURFACE,                            &
     &          NGNOD2D,                                                &
     &          STACEY_ABSORBING_CONDITIONS,                            &
     &          PML_CONDITIONS,                                         &
     &          NGLLX,                                                  &
     &          NGLLY,                                                  &
     &          NGLLZ,                                                  &
     &          NDIM,                                                   &
     &          NGNOD2D_FOUR_CORNERS,                                   &
     &          IMAIN,                                                  &
     &          BOTTOM_FREE_SURFACE
              INTEGER(KIND=4) :: NSPEC2D_TOP
              INTEGER(KIND=4) :: NSPEC2D_BOTTOM
              INTEGER(KIND=4) :: NSPEC2D_YMAX
              INTEGER(KIND=4) :: NSPEC2D_YMIN
              INTEGER(KIND=4) :: NSPEC2D_XMAX
              INTEGER(KIND=4) :: NSPEC2D_XMIN
              INTEGER(KIND=4) :: NNODES_EXT_MESH
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              REAL(KIND=8) :: NODES_COORDS_EXT_MESH(3,NNODES_EXT_MESH)
              INTEGER(KIND=4) :: IBELM_XMIN(NSPEC2D_XMIN)
              INTEGER(KIND=4) :: IBELM_XMAX(NSPEC2D_XMAX)
              INTEGER(KIND=4) :: IBELM_YMIN(NSPEC2D_YMIN)
              INTEGER(KIND=4) :: IBELM_YMAX(NSPEC2D_YMAX)
              INTEGER(KIND=4) :: IBELM_BOTTOM(NSPEC2D_BOTTOM)
              INTEGER(KIND=4) :: IBELM_TOP(NSPEC2D_TOP)
              INTEGER(KIND=4) :: NODES_IBELM_XMIN(NGNOD2D,NSPEC2D_XMIN)
              INTEGER(KIND=4) :: NODES_IBELM_XMAX(NGNOD2D,NSPEC2D_XMAX)
              INTEGER(KIND=4) :: NODES_IBELM_YMIN(NGNOD2D,NSPEC2D_YMIN)
              INTEGER(KIND=4) :: NODES_IBELM_YMAX(NGNOD2D,NSPEC2D_YMAX)
              INTEGER(KIND=4) :: NODES_IBELM_BOTTOM(NGNOD2D,            &
     &NSPEC2D_BOTTOM)
              INTEGER(KIND=4) :: NODES_IBELM_TOP(NGNOD2D,NSPEC2D_TOP)
            END SUBROUTINE GET_ABSORBING_BOUNDARY
          END INTERFACE 
        END MODULE GET_ABSORBING_BOUNDARY__genmod
