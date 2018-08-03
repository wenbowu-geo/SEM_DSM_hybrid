        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MULTIPLY_ARRAYS_SOURCE__genmod
          INTERFACE 
            SUBROUTINE MULTIPLY_ARRAYS_SOURCE(SOURCEARRAYD,G11,G12,G13, &
     &G21,G22,G23,G31,G32,G33,HXIS,HPXIS,HETAS,HPETAS,HGAMMAS,HPGAMMAS,K&
     &,L,M)
              REAL(KIND=8) :: SOURCEARRAYD(3,5,5,5)
              REAL(KIND=8) :: G11(5,5,5)
              REAL(KIND=8) :: G12(5,5,5)
              REAL(KIND=8) :: G13(5,5,5)
              REAL(KIND=8) :: G21(5,5,5)
              REAL(KIND=8) :: G22(5,5,5)
              REAL(KIND=8) :: G23(5,5,5)
              REAL(KIND=8) :: G31(5,5,5)
              REAL(KIND=8) :: G32(5,5,5)
              REAL(KIND=8) :: G33(5,5,5)
              REAL(KIND=8) :: HXIS(5)
              REAL(KIND=8) :: HPXIS(5)
              REAL(KIND=8) :: HETAS(5)
              REAL(KIND=8) :: HPETAS(5)
              REAL(KIND=8) :: HGAMMAS(5)
              REAL(KIND=8) :: HPGAMMAS(5)
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: L
              INTEGER(KIND=4) :: M
            END SUBROUTINE MULTIPLY_ARRAYS_SOURCE
          END INTERFACE 
        END MODULE MULTIPLY_ARRAYS_SOURCE__genmod
