        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:38 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_SU_HEADER__genmod
          INTERFACE 
            SUBROUTINE WRITE_SU_HEADER(IREC_LOCAL,IREC,NREC,NSTEP,DT,DX,&
     &X_FOUND,Y_FOUND,Z_FOUND,X_FOUND_SOURCE,Y_FOUND_SOURCE,            &
     &Z_FOUND_SOURCE)
              INTEGER(KIND=4) :: IREC_LOCAL
              INTEGER(KIND=4) :: IREC
              INTEGER(KIND=4) :: NREC
              INTEGER(KIND=4) :: NSTEP
              REAL(KIND=8) :: DT
              REAL(KIND=4) :: DX
              REAL(KIND=8) :: X_FOUND
              REAL(KIND=8) :: Y_FOUND
              REAL(KIND=8) :: Z_FOUND
              REAL(KIND=8) :: X_FOUND_SOURCE
              REAL(KIND=8) :: Y_FOUND_SOURCE
              REAL(KIND=8) :: Z_FOUND_SOURCE
            END SUBROUTINE WRITE_SU_HEADER
          END INTERFACE 
        END MODULE WRITE_SU_HEADER__genmod
