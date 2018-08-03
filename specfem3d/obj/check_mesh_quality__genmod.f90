        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:43 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_MESH_QUALITY__genmod
          INTERFACE 
            SUBROUTINE CHECK_MESH_QUALITY(MYRANK,VP_MAX,NGLOB,NSPEC,X,Y,&
     &Z,IBOOL,CREATE_VTK_FILES,PRNAME)
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: VP_MAX
              REAL(KIND=8) :: X(NGLOB)
              REAL(KIND=8) :: Y(NGLOB)
              REAL(KIND=8) :: Z(NGLOB)
              INTEGER(KIND=4) :: IBOOL(8,NSPEC)
              LOGICAL(KIND=4) :: CREATE_VTK_FILES
              CHARACTER(LEN=512) :: PRNAME
            END SUBROUTINE CHECK_MESH_QUALITY
          END INTERFACE 
        END MODULE CHECK_MESH_QUALITY__genmod
