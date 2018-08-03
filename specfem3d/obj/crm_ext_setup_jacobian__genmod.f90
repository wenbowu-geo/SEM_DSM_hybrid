        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:59 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CRM_EXT_SETUP_JACOBIAN__genmod
          INTERFACE 
            SUBROUTINE CRM_EXT_SETUP_JACOBIAN(MYRANK,XSTORE,YSTORE,     &
     &ZSTORE,NSPEC,NODES_COORDS_EXT_MESH,NNODES_EXT_MESH,ELMNTS_EXT_MESH&
     &,NELMNTS_EXT_MESH)
              USE GENERATE_DATABASES_PAR, ONLY :                        &
     &          NGNOD,                                                  &
     &          NGNOD2D,                                                &
     &          NDIM,                                                   &
     &          NGLLX,                                                  &
     &          NGLLY,                                                  &
     &          NGLLZ,                                                  &
     &          GAUSSALPHA,                                             &
     &          GAUSSBETA
              INTEGER(KIND=4) :: NELMNTS_EXT_MESH
              INTEGER(KIND=4) :: NNODES_EXT_MESH
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: XSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: YSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: ZSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: NODES_COORDS_EXT_MESH(3,NNODES_EXT_MESH)
              INTEGER(KIND=4) :: ELMNTS_EXT_MESH(NGNOD,NELMNTS_EXT_MESH)
            END SUBROUTINE CRM_EXT_SETUP_JACOBIAN
          END INTERFACE 
        END MODULE CRM_EXT_SETUP_JACOBIAN__genmod
