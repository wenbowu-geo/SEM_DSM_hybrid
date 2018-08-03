        !COMPILER-GENERATED INTERFACE MODULE: Fri Aug  3 10:43:32 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_ATTENUATION_PROPERTY_VALUES__genmod
          INTERFACE 
            SUBROUTINE GET_ATTENUATION_PROPERTY_VALUES(TAU_S,TAU_EPS,   &
     &BETA,ONE_MINUS_SUM_BETA)
              REAL(KIND=8), INTENT(IN) :: TAU_S(3)
              REAL(KIND=8), INTENT(IN) :: TAU_EPS(3)
              REAL(KIND=8), INTENT(OUT) :: BETA(3)
              REAL(KIND=8), INTENT(OUT) :: ONE_MINUS_SUM_BETA
            END SUBROUTINE GET_ATTENUATION_PROPERTY_VALUES
          END INTERFACE 
        END MODULE GET_ATTENUATION_PROPERTY_VALUES__genmod
