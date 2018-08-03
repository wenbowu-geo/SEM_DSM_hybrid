        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:53 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALC_JACOBIAN__genmod
          INTERFACE 
            SUBROUTINE CALC_JACOBIAN(MYRANK,XIXSTORE,XIYSTORE,XIZSTORE, &
     &ETAXSTORE,ETAYSTORE,ETAZSTORE,GAMMAXSTORE,GAMMAYSTORE,GAMMAZSTORE,&
     &JACOBIANSTORE,XSTORE,YSTORE,ZSTORE,XELM,YELM,ZELM,SHAPE3D,        &
     &DERSHAPE3D,ISPEC,NSPEC)
              USE GENERATE_DATABASES_PAR, ONLY :                        &
     &          NGNOD,                                                  &
     &          CUSTOM_REAL,                                            &
     &          NGLLX,                                                  &
     &          NGLLY,                                                  &
     &          NGLLZ,                                                  &
     &          NDIM,                                                   &
     &          SIZE_REAL,                                              &
     &          ZERO
              INTEGER(KIND=4) :: NSPEC
              INTEGER(KIND=4) :: MYRANK
              REAL(KIND=8) :: XIXSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: XIYSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: XIZSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: ETAXSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: ETAYSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: ETAZSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: GAMMAXSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: GAMMAYSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: GAMMAZSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: JACOBIANSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: XSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: YSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: ZSTORE(5,5,5,NSPEC)
              REAL(KIND=8) :: XELM(NGNOD)
              REAL(KIND=8) :: YELM(NGNOD)
              REAL(KIND=8) :: ZELM(NGNOD)
              REAL(KIND=8) :: SHAPE3D(NGNOD,5,5,5)
              REAL(KIND=8) :: DERSHAPE3D(3,NGNOD,5,5,5)
              INTEGER(KIND=4) :: ISPEC
            END SUBROUTINE CALC_JACOBIAN
          END INTERFACE 
        END MODULE CALC_JACOBIAN__genmod
