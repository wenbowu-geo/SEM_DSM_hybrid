        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:58 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_STRAIN_PRODUCT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_STRAIN_PRODUCT(PROD,EPS_TRACE_OVER_3,    &
     &EPSDEV,B_EPS_TRACE_OVER_3,B_EPSDEV)
              REAL(KIND=8) :: PROD(21)
              REAL(KIND=8) :: EPS_TRACE_OVER_3
              REAL(KIND=8) :: EPSDEV(5)
              REAL(KIND=8) :: B_EPS_TRACE_OVER_3
              REAL(KIND=8) :: B_EPSDEV(5)
            END SUBROUTINE COMPUTE_STRAIN_PRODUCT
          END INTERFACE 
        END MODULE COMPUTE_STRAIN_PRODUCT__genmod
