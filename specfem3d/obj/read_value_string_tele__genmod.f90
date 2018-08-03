        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:44:03 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_VALUE_STRING_TELE__genmod
          INTERFACE 
            SUBROUTINE READ_VALUE_STRING_TELE(IUNIT,IGNORE_JUNK,        &
     &VALUE_TO_READ,NAME,IER)
              INTEGER(KIND=4) :: IUNIT
              LOGICAL(KIND=4) :: IGNORE_JUNK
              CHARACTER(*) :: VALUE_TO_READ
              CHARACTER(*) :: NAME
              INTEGER(KIND=4) :: IER
            END SUBROUTINE READ_VALUE_STRING_TELE
          END INTERFACE 
        END MODULE READ_VALUE_STRING_TELE__genmod
