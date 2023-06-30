        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:01 2020
        MODULE CONSUP__genmod
          INTERFACE 
            SUBROUTINE CONSUP(UIN,FLUX,DIV,DT,NGRID)
              REAL(KIND=8) :: UIN(1:32,-1:4,-1:4,-1:4,1:8)
              REAL(KIND=8) :: FLUX(1:32,1:3,1:3,1:3,1:8,1:3)
              REAL(KIND=8) :: DIV(1:32,1:3,1:3,1:3)
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE CONSUP
          END INTERFACE 
        END MODULE CONSUP__genmod
