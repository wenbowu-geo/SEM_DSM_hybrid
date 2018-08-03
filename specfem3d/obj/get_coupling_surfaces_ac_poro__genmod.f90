        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_COUPLING_SURFACES_AC_PORO__genmod
          INTERFACE 
            SUBROUTINE GET_COUPLING_SURFACES_AC_PORO(MYRANK,NSPEC,IBOOL,&
     &ACOUSTIC_FLAG)
              USE CREATE_REGIONS_MESH_EXT_PAR
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              INTEGER(KIND=4) :: IBOOL(5,5,5,NSPEC)
              INTEGER(KIND=4) :: ACOUSTIC_FLAG(NGLOB_DUMMY)
            END SUBROUTINE GET_COUPLING_SURFACES_AC_PORO
          END INTERFACE 
        END MODULE GET_COUPLING_SURFACES_AC_PORO__genmod
