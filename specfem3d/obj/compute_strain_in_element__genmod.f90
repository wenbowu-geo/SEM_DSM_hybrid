        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:56 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_STRAIN_IN_ELEMENT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_STRAIN_IN_ELEMENT(TEMPX1_ATT,TEMPX2_ATT, &
     &TEMPX3_ATT,TEMPX1,TEMPX2,TEMPX3,TEMPY1_ATT,TEMPY2_ATT,TEMPY3_ATT, &
     &TEMPY1,TEMPY2,TEMPY3,TEMPZ1_ATT,TEMPZ2_ATT,TEMPZ3_ATT,TEMPZ1,     &
     &TEMPZ2,TEMPZ3,DUMMYX_LOC,DUMMYY_LOC,DUMMYZ_LOC,HPRIME_XX,HPRIME_YY&
     &,HPRIME_ZZ)
              REAL(KIND=8) :: TEMPX1_ATT(5,5,5)
              REAL(KIND=8) :: TEMPX2_ATT(5,5,5)
              REAL(KIND=8) :: TEMPX3_ATT(5,5,5)
              REAL(KIND=8) :: TEMPX1(5,5,5)
              REAL(KIND=8) :: TEMPX2(5,5,5)
              REAL(KIND=8) :: TEMPX3(5,5,5)
              REAL(KIND=8) :: TEMPY1_ATT(5,5,5)
              REAL(KIND=8) :: TEMPY2_ATT(5,5,5)
              REAL(KIND=8) :: TEMPY3_ATT(5,5,5)
              REAL(KIND=8) :: TEMPY1(5,5,5)
              REAL(KIND=8) :: TEMPY2(5,5,5)
              REAL(KIND=8) :: TEMPY3(5,5,5)
              REAL(KIND=8) :: TEMPZ1_ATT(5,5,5)
              REAL(KIND=8) :: TEMPZ2_ATT(5,5,5)
              REAL(KIND=8) :: TEMPZ3_ATT(5,5,5)
              REAL(KIND=8) :: TEMPZ1(5,5,5)
              REAL(KIND=8) :: TEMPZ2(5,5,5)
              REAL(KIND=8) :: TEMPZ3(5,5,5)
              REAL(KIND=8) :: DUMMYX_LOC(5,5,5)
              REAL(KIND=8) :: DUMMYY_LOC(5,5,5)
              REAL(KIND=8) :: DUMMYZ_LOC(5,5,5)
              REAL(KIND=8) :: HPRIME_XX(5,5)
              REAL(KIND=8) :: HPRIME_YY(5,5)
              REAL(KIND=8) :: HPRIME_ZZ(5,5)
            END SUBROUTINE COMPUTE_STRAIN_IN_ELEMENT
          END INTERFACE 
        END MODULE COMPUTE_STRAIN_IN_ELEMENT__genmod
