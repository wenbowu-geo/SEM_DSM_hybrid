        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:52 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_TOPO_BATHY_FILE__genmod
          INTERFACE 
            SUBROUTINE READ_TOPO_BATHY_FILE(ITOPO_BATHY,NX_TOPO,NY_TOPO)
              INTEGER(KIND=4) :: NY_TOPO
              INTEGER(KIND=4) :: NX_TOPO
              INTEGER(KIND=4) :: ITOPO_BATHY(NX_TOPO,NY_TOPO)
            END SUBROUTINE READ_TOPO_BATHY_FILE
          END INTERFACE 
        END MODULE READ_TOPO_BATHY_FILE__genmod
