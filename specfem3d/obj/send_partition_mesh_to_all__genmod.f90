        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:42 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SEND_PARTITION_MESH_TO_ALL__genmod
          INTERFACE 
            SUBROUTINE SEND_PARTITION_MESH_TO_ALL(MYRANK,IPART,NPART)
              USE MODULE_MESH
              INTEGER(KIND=4), INTENT(IN) :: MYRANK
              INTEGER(KIND=4), INTENT(IN) :: IPART(NSPEC_GLOB)
              INTEGER(KIND=4), INTENT(IN) :: NPART
            END SUBROUTINE SEND_PARTITION_MESH_TO_ALL
          END INTERFACE 
        END MODULE SEND_PARTITION_MESH_TO_ALL__genmod
