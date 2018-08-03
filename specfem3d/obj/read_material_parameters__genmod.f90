        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:49 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_MATERIAL_PARAMETERS__genmod
          INTERFACE 
            SUBROUTINE READ_MATERIAL_PARAMETERS(IUNIT,MAT_ID,RHO,VP,VS, &
     &Q_FLAG,ANISOTROPY_FLAG,DOMAIN_ID,IER)
              INTEGER(KIND=4) :: IUNIT
              INTEGER(KIND=4) :: MAT_ID
              REAL(KIND=8) :: RHO
              REAL(KIND=8) :: VP
              REAL(KIND=8) :: VS
              REAL(KIND=8) :: Q_FLAG
              REAL(KIND=8) :: ANISOTROPY_FLAG
              INTEGER(KIND=4) :: DOMAIN_ID
              INTEGER(KIND=4) :: IER
            END SUBROUTINE READ_MATERIAL_PARAMETERS
          END INTERFACE 
        END MODULE READ_MATERIAL_PARAMETERS__genmod
