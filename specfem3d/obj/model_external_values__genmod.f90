        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:19 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MODEL_EXTERNAL_VALUES__genmod
          INTERFACE 
            SUBROUTINE MODEL_EXTERNAL_VALUES(XMESH,YMESH,ZMESH,RHO,VP,VS&
     &,QKAPPA_ATTEN,QMU_ATTEN,IFLAG_ANISO,IDOMAIN_ID)
              REAL(KIND=8), INTENT(IN) :: XMESH
              REAL(KIND=8), INTENT(IN) :: YMESH
              REAL(KIND=8), INTENT(IN) :: ZMESH
              REAL(KIND=8) :: RHO
              REAL(KIND=8) :: VP
              REAL(KIND=8) :: VS
              REAL(KIND=8) :: QKAPPA_ATTEN
              REAL(KIND=8) :: QMU_ATTEN
              INTEGER(KIND=4) :: IFLAG_ANISO
              INTEGER(KIND=4) :: IDOMAIN_ID
            END SUBROUTINE MODEL_EXTERNAL_VALUES
          END INTERFACE 
        END MODULE MODEL_EXTERNAL_VALUES__genmod
