        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:49 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CREATE_MESH_QUALITY_DATA_3D__genmod
          INTERFACE 
            SUBROUTINE CREATE_MESH_QUALITY_DATA_3D(X,Y,Z,IBOOL,ISPEC,   &
     &NSPEC,NGLOB,VP_MAX,DELTA_T,EQUIANGLE_SKEWNESS,EDGE_ASPECT_RATIO,  &
     &DIAGONAL_ASPECT_RATIO,STABILITY,DISTMIN,DISTMAX,DISTMEAN)
              INTEGER(KIND=4) :: NGLOB
              INTEGER(KIND=4) :: NSPEC
              REAL(KIND=8) :: X(NGLOB)
              REAL(KIND=8) :: Y(NGLOB)
              REAL(KIND=8) :: Z(NGLOB)
              INTEGER(KIND=4) :: IBOOL(8,NSPEC)
              INTEGER(KIND=4) :: ISPEC
              REAL(KIND=8) :: VP_MAX
              REAL(KIND=8) :: DELTA_T
              REAL(KIND=8) :: EQUIANGLE_SKEWNESS
              REAL(KIND=8) :: EDGE_ASPECT_RATIO
              REAL(KIND=8) :: DIAGONAL_ASPECT_RATIO
              REAL(KIND=8) :: STABILITY
              REAL(KIND=8) :: DISTMIN
              REAL(KIND=8) :: DISTMAX
              REAL(KIND=8) :: DISTMEAN
            END SUBROUTINE CREATE_MESH_QUALITY_DATA_3D
          END INTERFACE 
        END MODULE CREATE_MESH_QUALITY_DATA_3D__genmod
