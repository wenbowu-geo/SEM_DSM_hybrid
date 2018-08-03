        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_COUPLING_SURFACES__genmod
          INTERFACE 
            SUBROUTINE GET_COUPLING_SURFACES(MYRANK,NSPEC,IBOOL,NPROC,  &
     &NIBOOL_INTERFACES_EXT_MESH,IBOOL_INTERFACES_EXT_MESH,             &
     &NUM_INTERFACES_EXT_MESH,MAX_INTERFACE_SIZE_EXT_MESH,              &
     &MY_NEIGHBOURS_EXT_MESH)
              INTEGER(KIND=4) :: MAX_INTERFACE_SIZE_EXT_MESH
              INTEGER(KIND=4) :: NUM_INTERFACES_EXT_MESH
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              INTEGER(KIND=4) :: NPROC
              INTEGER(KIND=4) :: NIBOOL_INTERFACES_EXT_MESH(            &
     &NUM_INTERFACES_EXT_MESH)
              INTEGER(KIND=4) :: IBOOL_INTERFACES_EXT_MESH(5*5*         &
     &MAX_INTERFACE_SIZE_EXT_MESH,NUM_INTERFACES_EXT_MESH)
              INTEGER(KIND=4) :: MY_NEIGHBOURS_EXT_MESH(                &
     &NUM_INTERFACES_EXT_MESH)
            END SUBROUTINE GET_COUPLING_SURFACES
          END INTERFACE 
        END MODULE GET_COUPLING_SURFACES__genmod
