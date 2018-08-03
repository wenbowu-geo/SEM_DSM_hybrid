        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:30 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DEFINE_DERIVATION_MATRICES__genmod
          INTERFACE 
            SUBROUTINE DEFINE_DERIVATION_MATRICES(XIGLL,YIGLL,ZIGLL,    &
     &WXGLL,WYGLL,WZGLL,HPRIME_XX,HPRIME_YY,HPRIME_ZZ,HPRIMEWGLL_XX,    &
     &HPRIMEWGLL_YY,HPRIMEWGLL_ZZ,WGLLWGLL_XY,WGLLWGLL_XZ,WGLLWGLL_YZ)
              REAL(KIND=8) :: XIGLL(5)
              REAL(KIND=8) :: YIGLL(5)
              REAL(KIND=8) :: ZIGLL(5)
              REAL(KIND=8) :: WXGLL(5)
              REAL(KIND=8) :: WYGLL(5)
              REAL(KIND=8) :: WZGLL(5)
              REAL(KIND=8) :: HPRIME_XX(5,5)
              REAL(KIND=8) :: HPRIME_YY(5,5)
              REAL(KIND=8) :: HPRIME_ZZ(5,5)
              REAL(KIND=8) :: HPRIMEWGLL_XX(5,5)
              REAL(KIND=8) :: HPRIMEWGLL_YY(5,5)
              REAL(KIND=8) :: HPRIMEWGLL_ZZ(5,5)
              REAL(KIND=8) :: WGLLWGLL_XY(5,5)
              REAL(KIND=8) :: WGLLWGLL_XZ(5,5)
              REAL(KIND=8) :: WGLLWGLL_YZ(5,5)
            END SUBROUTINE DEFINE_DERIVATION_MATRICES
          END INTERFACE 
        END MODULE DEFINE_DERIVATION_MATRICES__genmod
