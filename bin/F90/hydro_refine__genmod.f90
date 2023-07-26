        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:02 2020
        MODULE HYDRO_REFINE__genmod
          INTERFACE 
            SUBROUTINE HYDRO_REFINE(UG,UM,UD,OK,NN)
              REAL(KIND=8) :: UG(1:32,1:8)
              REAL(KIND=8) :: UM(1:32,1:8)
              REAL(KIND=8) :: UD(1:32,1:8)
              LOGICAL(KIND=4) :: OK(1:32)
              INTEGER(KIND=4) :: NN
            END SUBROUTINE HYDRO_REFINE
          END INTERFACE 
        END MODULE HYDRO_REFINE__genmod
