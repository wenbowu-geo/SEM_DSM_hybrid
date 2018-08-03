        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MODEL_DEFAULT__genmod
          INTERFACE 
            SUBROUTINE MODEL_DEFAULT(MATERIALS_EXT_MESH,NMAT_EXT_MESH,  &
     &UNDEF_MAT_PROP,NUNDEFMAT_EXT_MESH,IMATERIAL_ID,IMATERIAL_DEF,XMESH&
     &,YMESH,ZMESH,RHO,VP,VS,IFLAG_ANISO,QKAPPA_ATTEN,QMU_ATTEN,        &
     &IDOMAIN_ID,RHO_S,KAPPA_S,RHO_F,KAPPA_F,ETA_F,KAPPA_FR,MU_FR,PHI,  &
     &TORT,KXX,KXY,KXZ,KYY,KYZ,KZZ)
              INTEGER(KIND=4), INTENT(IN) :: NUNDEFMAT_EXT_MESH
              INTEGER(KIND=4), INTENT(IN) :: NMAT_EXT_MESH
              REAL(KIND=8), INTENT(IN) :: MATERIALS_EXT_MESH(16,        &
     &NMAT_EXT_MESH)
              CHARACTER(LEN=512) :: UNDEF_MAT_PROP(6,NUNDEFMAT_EXT_MESH)
              INTEGER(KIND=4), INTENT(IN) :: IMATERIAL_ID
              INTEGER(KIND=4), INTENT(IN) :: IMATERIAL_DEF
              REAL(KIND=8), INTENT(IN) :: XMESH
              REAL(KIND=8), INTENT(IN) :: YMESH
              REAL(KIND=8), INTENT(IN) :: ZMESH
              REAL(KIND=8) :: RHO
              REAL(KIND=8) :: VP
              REAL(KIND=8) :: VS
              INTEGER(KIND=4) :: IFLAG_ANISO
              REAL(KIND=8) :: QKAPPA_ATTEN
              REAL(KIND=8) :: QMU_ATTEN
              INTEGER(KIND=4) :: IDOMAIN_ID
              REAL(KIND=8) :: RHO_S
              REAL(KIND=8) :: KAPPA_S
              REAL(KIND=8) :: RHO_F
              REAL(KIND=8) :: KAPPA_F
              REAL(KIND=8) :: ETA_F
              REAL(KIND=8) :: KAPPA_FR
              REAL(KIND=8) :: MU_FR
              REAL(KIND=8) :: PHI
              REAL(KIND=8) :: TORT
              REAL(KIND=8) :: KXX
              REAL(KIND=8) :: KXY
              REAL(KIND=8) :: KXZ
              REAL(KIND=8) :: KYY
              REAL(KIND=8) :: KYZ
              REAL(KIND=8) :: KZZ
            END SUBROUTINE MODEL_DEFAULT
          END INTERFACE 
        END MODULE MODEL_DEFAULT__genmod
