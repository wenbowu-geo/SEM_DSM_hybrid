        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:17 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_MODEL_VALUES__genmod
          INTERFACE 
            SUBROUTINE GET_MODEL_VALUES(MATERIALS_EXT_MESH,NMAT_EXT_MESH&
     &,UNDEF_MAT_PROP,NUNDEFMAT_EXT_MESH,IMATERIAL_ID,IMATERIAL_DEF,    &
     &XMESH,YMESH,ZMESH,XMESH_CUBEDSPH,YMESH_CUBEDSPH,ZMESH_CUBEDSPH,   &
     &R_MIDDLE,RHO,VP,VS,QKAPPA_ATTEN,QMU_ATTEN,IDOMAIN_ID,RHO_S,KAPPA_S&
     &,RHO_F,KAPPA_F,ETA_F,KAPPA_FR,MU_FR,PHI,TORT,KXX,KXY,KXZ,KYY,KYZ, &
     &KZZ,C11,C12,C13,C14,C15,C16,C22,C23,C24,C25,C26,C33,C34,C35,C36,  &
     &C44,C45,C46,C55,C56,C66,ANISOTROPY,MYRANK)
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
              REAL(KIND=8), INTENT(IN) :: XMESH_CUBEDSPH
              REAL(KIND=8), INTENT(IN) :: YMESH_CUBEDSPH
              REAL(KIND=8), INTENT(IN) :: ZMESH_CUBEDSPH
              REAL(KIND=8), INTENT(IN) :: R_MIDDLE
              REAL(KIND=8) :: RHO
              REAL(KIND=8) :: VP
              REAL(KIND=8) :: VS
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
              REAL(KIND=8) :: C11
              REAL(KIND=8) :: C12
              REAL(KIND=8) :: C13
              REAL(KIND=8) :: C14
              REAL(KIND=8) :: C15
              REAL(KIND=8) :: C16
              REAL(KIND=8) :: C22
              REAL(KIND=8) :: C23
              REAL(KIND=8) :: C24
              REAL(KIND=8) :: C25
              REAL(KIND=8) :: C26
              REAL(KIND=8) :: C33
              REAL(KIND=8) :: C34
              REAL(KIND=8) :: C35
              REAL(KIND=8) :: C36
              REAL(KIND=8) :: C44
              REAL(KIND=8) :: C45
              REAL(KIND=8) :: C46
              REAL(KIND=8) :: C55
              REAL(KIND=8) :: C56
              REAL(KIND=8) :: C66
              LOGICAL(KIND=4) :: ANISOTROPY
              INTEGER(KIND=4) :: MYRANK
            END SUBROUTINE GET_MODEL_VALUES
          END INTERFACE 
        END MODULE GET_MODEL_VALUES__genmod
