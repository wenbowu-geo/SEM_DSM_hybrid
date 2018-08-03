        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:17 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_MPI__genmod
          INTERFACE 
            SUBROUTINE GET_MPI(MYRANK,NGLOB,NSPEC,IBOOL,NELMNTS_EXT_MESH&
     &,ELMNTS_EXT_MESH,MY_NELMNTS_NEIGHBOURS_EXT_MESH,                  &
     &MY_INTERFACES_EXT_MESH,IBOOL_INTERFACES_EXT_MESH,                 &
     &NIBOOL_INTERFACES_EXT_MESH,NUM_INTERFACES_EXT_MESH,               &
     &MAX_INTERFACE_SIZE_EXT_MESH,MY_NEIGHBOURS_EXT_MESH)
              USE GENERATE_DATABASES_PAR, ONLY :                        &
     &          NPROC,                                                  &
     &          NGNOD,                                                  &
     &          NGLLX,                                                  &
     &          NGLLY,                                                  &
     &          NGLLZ,                                                  &
     &          SMALLVAL_TOL,                                           &
     &          IMAIN
              INTEGER(KIND=4) :: MAX_INTERFACE_SIZE_EXT_MESH
              INTEGER(KIND=4) :: NUM_INTERFACES_EXT_MESH
              INTEGER(KIND=4) :: NELMNTS_EXT_MESH
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              INTEGER(KIND=4) :: ELMNTS_EXT_MESH(NGNOD,NELMNTS_EXT_MESH)
              INTEGER(KIND=4) :: MY_NELMNTS_NEIGHBOURS_EXT_MESH(        &
     &NUM_INTERFACES_EXT_MESH)
              INTEGER(KIND=4) :: MY_INTERFACES_EXT_MESH(6,              &
     &MAX_INTERFACE_SIZE_EXT_MESH,NUM_INTERFACES_EXT_MESH)
              INTEGER(KIND=4) :: IBOOL_INTERFACES_EXT_MESH(5*5*         &
     &MAX_INTERFACE_SIZE_EXT_MESH,NUM_INTERFACES_EXT_MESH)
              INTEGER(KIND=4) :: NIBOOL_INTERFACES_EXT_MESH(            &
     &NUM_INTERFACES_EXT_MESH)
              INTEGER(KIND=4) :: MY_NEIGHBOURS_EXT_MESH(                &
     &NUM_INTERFACES_EXT_MESH)
            END SUBROUTINE GET_MPI
          END INTERFACE 
        END MODULE GET_MPI__genmod
