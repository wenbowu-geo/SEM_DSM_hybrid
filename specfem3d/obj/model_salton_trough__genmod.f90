        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:20 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MODEL_SALTON_TROUGH__genmod
          INTERFACE 
            SUBROUTINE MODEL_SALTON_TROUGH(XMESH,YMESH,ZMESH,RHO,VP,VS, &
     &QMU_ATTEN)
              REAL(KIND=8), INTENT(IN) :: XMESH
              REAL(KIND=8), INTENT(IN) :: YMESH
              REAL(KIND=8), INTENT(IN) :: ZMESH
              REAL(KIND=8) :: RHO
              REAL(KIND=8) :: VP
              REAL(KIND=8) :: VS
              REAL(KIND=8) :: QMU_ATTEN
            END SUBROUTINE MODEL_SALTON_TROUGH
          END INTERFACE 
        END MODULE MODEL_SALTON_TROUGH__genmod
