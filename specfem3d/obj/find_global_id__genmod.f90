        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:42:57 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FIND_GLOBAL_ID__genmod
          INTERFACE 
            SUBROUTINE FIND_GLOBAL_ID(INPUT_ARRAY,SIZE_ARRAY,           &
     &NPOINTS_ARRAY,TOLERENCE_VALUE,ID_INGLOBAL,TABLE_GLOBAL,           &
     &N_TABLE_GLOBAL,MAX_NTABLE_GLOBAL,NPROC,MYRANK)
              INTEGER(KIND=4), INTENT(IN) :: MAX_NTABLE_GLOBAL
              INTEGER(KIND=4), INTENT(IN) :: SIZE_ARRAY
              REAL(KIND=8), INTENT(IN) :: INPUT_ARRAY(SIZE_ARRAY)
              INTEGER(KIND=4), INTENT(IN) :: NPOINTS_ARRAY
              REAL(KIND=8), INTENT(IN) :: TOLERENCE_VALUE
              INTEGER(KIND=4), INTENT(OUT) :: ID_INGLOBAL(SIZE_ARRAY)
              REAL(KIND=8), INTENT(OUT) :: TABLE_GLOBAL(                &
     &MAX_NTABLE_GLOBAL)
              INTEGER(KIND=4), INTENT(OUT) :: N_TABLE_GLOBAL
              INTEGER(KIND=4), INTENT(IN) :: NPROC
              INTEGER(KIND=4), INTENT(IN) :: MYRANK
            END SUBROUTINE FIND_GLOBAL_ID
          END INTERFACE 
        END MODULE FIND_GLOBAL_ID__genmod
