        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:12 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DERIV__genmod
          INTERFACE 
            SUBROUTINE DERIV(Y,YPRIME,N,R,NDIS,KDIS,S1,S2,S3)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: YPRIME(N)
              REAL(KIND=8) :: R(N)
              INTEGER(KIND=4) :: NDIS
              INTEGER(KIND=4) :: KDIS(28)
              REAL(KIND=8) :: S1(N)
              REAL(KIND=8) :: S2(N)
              REAL(KIND=8) :: S3(N)
            END SUBROUTINE DERIV
          END INTERFACE 
        END MODULE DERIV__genmod
