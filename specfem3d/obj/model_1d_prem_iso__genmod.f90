        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MODEL_1D_PREM_ISO__genmod
          INTERFACE 
            SUBROUTINE MODEL_1D_PREM_ISO(XMESH,YMESH,ZMESH,RHO_PREM,    &
     &VP_PREM,VS_PREM,QMU_ATTEN)
              REAL(KIND=8), INTENT(IN) :: XMESH
              REAL(KIND=8), INTENT(IN) :: YMESH
              REAL(KIND=8), INTENT(IN) :: ZMESH
              REAL(KIND=8), INTENT(INOUT) :: RHO_PREM
              REAL(KIND=8), INTENT(INOUT) :: VP_PREM
              REAL(KIND=8), INTENT(INOUT) :: VS_PREM
              REAL(KIND=8), INTENT(INOUT) :: QMU_ATTEN
            END SUBROUTINE MODEL_1D_PREM_ISO
          END INTERFACE 
        END MODULE MODEL_1D_PREM_ISO__genmod
