        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:45 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CVD_COUNT_TOTALS_EXT_MESH__genmod
          INTERFACE 
            SUBROUTINE CVD_COUNT_TOTALS_EXT_MESH(NUM_NODE,NODE_LIST,    &
     &LOCAL_PATH,NPP,NEE,HIGH_RESOLUTION_MESH,MESH_HANDLE,ADIOS_FOR_MESH&
     &)
              INTEGER(KIND=4), INTENT(IN) :: NUM_NODE
              INTEGER(KIND=4), INTENT(IN) :: NODE_LIST(600)
              CHARACTER(LEN=512), INTENT(IN) :: LOCAL_PATH
              INTEGER(KIND=4), INTENT(OUT) :: NPP
              INTEGER(KIND=4), INTENT(OUT) :: NEE
              LOGICAL(KIND=4), INTENT(IN) :: HIGH_RESOLUTION_MESH
              INTEGER(KIND=8), INTENT(IN) :: MESH_HANDLE
              LOGICAL(KIND=4), INTENT(IN) :: ADIOS_FOR_MESH
            END SUBROUTINE CVD_COUNT_TOTALS_EXT_MESH
          END INTERFACE 
        END MODULE CVD_COUNT_TOTALS_EXT_MESH__genmod
