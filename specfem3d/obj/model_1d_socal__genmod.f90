        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MODEL_1D_SOCAL__genmod
          INTERFACE 
            SUBROUTINE MODEL_1D_SOCAL(XMESH,YMESH,ZMESH,RHO,VP,VS,      &
     &QMU_ATTEN)
              REAL(KIND=8), INTENT(IN) :: XMESH
              REAL(KIND=8), INTENT(IN) :: YMESH
              REAL(KIND=8), INTENT(IN) :: ZMESH
              REAL(KIND=8), INTENT(INOUT) :: RHO
              REAL(KIND=8), INTENT(INOUT) :: VP
              REAL(KIND=8), INTENT(INOUT) :: VS
              REAL(KIND=8), INTENT(INOUT) :: QMU_ATTEN
            END SUBROUTINE MODEL_1D_SOCAL
          END INTERFACE 
        END MODULE MODEL_1D_SOCAL__genmod
