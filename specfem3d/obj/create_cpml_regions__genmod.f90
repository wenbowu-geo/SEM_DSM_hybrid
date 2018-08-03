        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CREATE_CPML_REGIONS__genmod
          INTERFACE 
            SUBROUTINE CREATE_CPML_REGIONS(NSPEC,NGLOB,NODES_COORDS)
              INTEGER(KIND=4), INTENT(IN) :: NGLOB
              INTEGER(KIND=4), INTENT(IN) :: NSPEC
              REAL(KIND=8), INTENT(IN) :: NODES_COORDS(NGLOB,3)
            END SUBROUTINE CREATE_CPML_REGIONS
          END INTERFACE 
        END MODULE CREATE_CPML_REGIONS__genmod
