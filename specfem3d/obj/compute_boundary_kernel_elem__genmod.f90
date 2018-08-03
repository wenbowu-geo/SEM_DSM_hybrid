        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:42 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_BOUNDARY_KERNEL_ELEM__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_BOUNDARY_KERNEL_ELEM(KERNEL,MUL,KAPPAL,  &
     &RHO_VSL,ACCEL,B_DISPL,DS,B_DS,NORM)
              REAL(KIND=8) :: KERNEL
              REAL(KIND=8) :: MUL
              REAL(KIND=8) :: KAPPAL
              REAL(KIND=8) :: RHO_VSL
              REAL(KIND=8) :: ACCEL(3)
              REAL(KIND=8) :: B_DISPL(3)
              REAL(KIND=8) :: DS(3,3)
              REAL(KIND=8) :: B_DS(3,3)
              REAL(KIND=8) :: NORM(3)
            END SUBROUTINE COMPUTE_BOUNDARY_KERNEL_ELEM
          END INTERFACE 
        END MODULE COMPUTE_BOUNDARY_KERNEL_ELEM__genmod
