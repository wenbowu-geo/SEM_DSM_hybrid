        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:56 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ADD_SLAB__genmod
          INTERFACE 
            SUBROUTINE ADD_SLAB(X,Y,Z,X_CUBEDSPH,Y_CUBEDSPH,Z_CUBEDSPH, &
     &VP,VS,RHO,MYRANK)
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y
              REAL(KIND=8) :: Z
              REAL(KIND=8) :: X_CUBEDSPH
              REAL(KIND=8) :: Y_CUBEDSPH
              REAL(KIND=8) :: Z_CUBEDSPH
              REAL(KIND=8) :: VP
              REAL(KIND=8) :: VS
              REAL(KIND=8) :: RHO
              INTEGER(KIND=4) :: MYRANK
            END SUBROUTINE ADD_SLAB
          END INTERFACE 
        END MODULE ADD_SLAB__genmod
