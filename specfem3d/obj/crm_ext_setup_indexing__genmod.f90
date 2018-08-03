        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:59 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CRM_EXT_SETUP_INDEXING__genmod
          INTERFACE 
            SUBROUTINE CRM_EXT_SETUP_INDEXING(IBOOL,XSTORE,YSTORE,ZSTORE&
     &,NSPEC,NGLOB,NPOINTOT,NNODES_EXT_MESH,NODES_COORDS_EXT_MESH,MYRANK&
     &)
              INTEGER(KIND=4) :: NNODES_EXT_MESH
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              REAL(KIND=8) :: XSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: YSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: ZSTORE(5,5,5,NSPEC)
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: NPOINTOT
              REAL(KIND=8) :: NODES_COORDS_EXT_MESH(3,NNODES_EXT_MESH)
              INTEGER(KIND=4) :: MYRANK
            END SUBROUTINE CRM_EXT_SETUP_INDEXING
          END INTERFACE 
        END MODULE CRM_EXT_SETUP_INDEXING__genmod
