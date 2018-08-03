        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:01 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_PNM_INITIALIZE__genmod
          INTERFACE 
            SUBROUTINE WRITE_PNM_INITIALIZE
              USE SPECFEM_PAR, ONLY :                                   &
     &          NGLOB_AB,                                               &
     &          NSPEC_AB,                                               &
     &          NPROC,                                                  &
     &          IBOOL,                                                  &
     &          XSTORE,                                                 &
     &          YSTORE,                                                 &
     &          ZSTORE,                                                 &
     &          NUM_INTERFACES_EXT_MESH,                                &
     &          MAX_NIBOOL_INTERFACES_EXT_MESH,                         &
     &          NIBOOL_INTERFACES_EXT_MESH,                             &
     &          MY_NEIGHBOURS_EXT_MESH,                                 &
     &          IBOOL_INTERFACES_EXT_MESH,                              &
     &          MYRANK
            END SUBROUTINE WRITE_PNM_INITIALIZE
          END INTERFACE 
        END MODULE WRITE_PNM_INITIALIZE__genmod
