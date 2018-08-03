        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:29 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SAVE_KERNELS_ELASTIC__genmod
          INTERFACE 
            SUBROUTINE SAVE_KERNELS_ELASTIC(ADIOS_HANDLE,ALPHAV_KL,     &
     &ALPHAH_KL,BETAV_KL,BETAH_KL,ETA_KL,RHOP_KL,ALPHA_KL,BETA_KL)
              INTEGER(KIND=8) :: ADIOS_HANDLE
              REAL(KIND=8) ,ALLOCATABLE :: ALPHAV_KL(:,:,:,:)
              REAL(KIND=8) ,ALLOCATABLE :: ALPHAH_KL(:,:,:,:)
              REAL(KIND=8) ,ALLOCATABLE :: BETAV_KL(:,:,:,:)
              REAL(KIND=8) ,ALLOCATABLE :: BETAH_KL(:,:,:,:)
              REAL(KIND=8) ,ALLOCATABLE :: ETA_KL(:,:,:,:)
              REAL(KIND=8) ,ALLOCATABLE :: RHOP_KL(:,:,:,:)
              REAL(KIND=8) ,ALLOCATABLE :: ALPHA_KL(:,:,:,:)
              REAL(KIND=8) ,ALLOCATABLE :: BETA_KL(:,:,:,:)
            END SUBROUTINE SAVE_KERNELS_ELASTIC
          END INTERFACE 
        END MODULE SAVE_KERNELS_ELASTIC__genmod
