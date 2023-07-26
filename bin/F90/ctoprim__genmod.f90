        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:01 2020
        MODULE CTOPRIM__genmod
          INTERFACE 
            SUBROUTINE CTOPRIM(UIN,Q,C,GRAVIN,DT,NGRID)
              REAL(KIND=8) :: UIN(1:32,-1:4,-1:4,-1:4,1:8)
              REAL(KIND=8) :: Q(1:32,-1:4,-1:4,-1:4,1:8)
              REAL(KIND=8) :: C(1:32,-1:4,-1:4,-1:4)
              REAL(KIND=8) :: GRAVIN(1:32,-1:4,-1:4,-1:4,1:3)
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE CTOPRIM
          END INTERFACE 
        END MODULE CTOPRIM__genmod
