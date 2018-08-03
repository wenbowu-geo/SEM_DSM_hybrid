        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:59 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CRM_SETUP_MOHO__genmod
          INTERFACE 
            SUBROUTINE CRM_SETUP_MOHO(MYRANK,NSPEC,NSPEC2D_MOHO_EXT,    &
     &IBELM_MOHO,NODES_IBELM_MOHO,NODES_COORDS_EXT_MESH,NNODES_EXT_MESH,&
     &IBOOL)
              USE GENERATE_DATABASES_PAR, ONLY :                        &
     &          NGNOD2D,                                                &
     &          NGLLX,                                                  &
     &          NGLLY,                                                  &
     &          NGLLZ,                                                  &
     &          CUSTOM_REAL,                                            &
     &          SIZE_REAL,                                              &
     &          IMAIN,                                                  &
     &          NDIM,                                                   &
     &          NGLLSQUARE,                                             &
     &          NGNOD2D_FOUR_CORNERS
              INTEGER(KIND=4) :: NNODES_EXT_MESH
              INTEGER(KIND=4) :: NSPEC2D_MOHO_EXT
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: IBELM_MOHO(NSPEC2D_MOHO_EXT)
              INTEGER(KIND=4) :: NODES_IBELM_MOHO(NGNOD2D,              &
     &NSPEC2D_MOHO_EXT)
              REAL(KIND=8) :: NODES_COORDS_EXT_MESH(3,NNODES_EXT_MESH)
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
            END SUBROUTINE CRM_SETUP_MOHO
          END INTERFACE 
        END MODULE CRM_SETUP_MOHO__genmod
