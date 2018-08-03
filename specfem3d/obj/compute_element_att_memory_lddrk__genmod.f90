        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:57 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_ELEMENT_ATT_MEMORY_LDDRK__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_ELEMENT_ATT_MEMORY_LDDRK(ISPEC,DELTAT,   &
     &NSPEC_AB,KAPPASTORE,MUSTORE,NSPEC_ATTENUATION_AB,                 &
     &FACTOR_COMMON_KAPPA,R_TRACE,EPSILONDEV_TRACE_LOC,                 &
     &NSPEC_ATTENUATION_AB_LDDRK,R_TRACE_LDDRK,FACTOR_COMMON,R_XX,R_YY, &
     &R_XY,R_XZ,R_YZ,R_XX_LDDRK,R_YY_LDDRK,R_XY_LDDRK,R_XZ_LDDRK,       &
     &R_YZ_LDDRK,EPSILONDEV_XX_LOC,EPSILONDEV_YY_LOC,EPSILONDEV_XY_LOC, &
     &EPSILONDEV_XZ_LOC,EPSILONDEV_YZ_LOC)
              INTEGER(KIND=4) :: NSPEC_ATTENUATION_AB_LDDRK
              INTEGER(KIND=4) :: NSPEC_ATTENUATION_AB
              INTEGER(KIND=4) :: NSPEC_AB
              INTEGER(KIND=4) :: ISPEC
              REAL(KIND=8) :: DELTAT
              REAL(KIND=8) :: KAPPASTORE(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: MUSTORE(5,5,5,NSPEC_AB)
              REAL(KIND=8) :: FACTOR_COMMON_KAPPA(3,5,5,5,              &
     &NSPEC_ATTENUATION_AB)
              REAL(KIND=8) :: R_TRACE(5,5,5,NSPEC_ATTENUATION_AB,3)
              REAL(KIND=8) :: EPSILONDEV_TRACE_LOC(5,5,5)
              REAL(KIND=8) :: R_TRACE_LDDRK(5,5,5,                      &
     &NSPEC_ATTENUATION_AB_LDDRK,3)
              REAL(KIND=8) :: FACTOR_COMMON(3,5,5,5,NSPEC_ATTENUATION_AB&
     &)
              REAL(KIND=8) :: R_XX(5,5,5,NSPEC_ATTENUATION_AB,3)
              REAL(KIND=8) :: R_YY(5,5,5,NSPEC_ATTENUATION_AB,3)
              REAL(KIND=8) :: R_XY(5,5,5,NSPEC_ATTENUATION_AB,3)
              REAL(KIND=8) :: R_XZ(5,5,5,NSPEC_ATTENUATION_AB,3)
              REAL(KIND=8) :: R_YZ(5,5,5,NSPEC_ATTENUATION_AB,3)
              REAL(KIND=8) :: R_XX_LDDRK(5,5,5,                         &
     &NSPEC_ATTENUATION_AB_LDDRK,3)
              REAL(KIND=8) :: R_YY_LDDRK(5,5,5,                         &
     &NSPEC_ATTENUATION_AB_LDDRK,3)
              REAL(KIND=8) :: R_XY_LDDRK(5,5,5,                         &
     &NSPEC_ATTENUATION_AB_LDDRK,3)
              REAL(KIND=8) :: R_XZ_LDDRK(5,5,5,                         &
     &NSPEC_ATTENUATION_AB_LDDRK,3)
              REAL(KIND=8) :: R_YZ_LDDRK(5,5,5,                         &
     &NSPEC_ATTENUATION_AB_LDDRK,3)
              REAL(KIND=8) :: EPSILONDEV_XX_LOC(5,5,5)
              REAL(KIND=8) :: EPSILONDEV_YY_LOC(5,5,5)
              REAL(KIND=8) :: EPSILONDEV_XY_LOC(5,5,5)
              REAL(KIND=8) :: EPSILONDEV_XZ_LOC(5,5,5)
              REAL(KIND=8) :: EPSILONDEV_YZ_LOC(5,5,5)
            END SUBROUTINE COMPUTE_ELEMENT_ATT_MEMORY_LDDRK
          END INTERFACE 
        END MODULE COMPUTE_ELEMENT_ATT_MEMORY_LDDRK__genmod
