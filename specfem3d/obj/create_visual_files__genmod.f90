        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:47 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CREATE_VISUAL_FILES__genmod
          INTERFACE 
            SUBROUTINE CREATE_VISUAL_FILES(CREATE_ABAQUS_FILES,         &
     &CREATE_DX_FILES,CREATE_VTK_FILES,NSPEC,NGLOB,PRNAME,NODES_COORDS, &
     &IBOOL,ISPEC_MATERIAL_ID)
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: NSPEC
              LOGICAL(KIND=4) :: CREATE_ABAQUS_FILES
              LOGICAL(KIND=4) :: CREATE_DX_FILES
              LOGICAL(KIND=4) :: CREATE_VTK_FILES
              CHARACTER(LEN=512) :: PRNAME
              REAL(KIND=8) :: NODES_COORDS(NGLOB,3)
              INTEGER(KIND=4) :: IBOOL(2,2,2,NSPEC)
              INTEGER(KIND=4) :: ISPEC_MATERIAL_ID(NSPEC)
            END SUBROUTINE CREATE_VISUAL_FILES
          END INTERFACE 
        END MODULE CREATE_VISUAL_FILES__genmod
