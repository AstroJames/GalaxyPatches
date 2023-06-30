        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:01 2020
        MODULE USLOPE__genmod
          INTERFACE 
            SUBROUTINE USLOPE(Q,DQ,DX,DT,NGRID)
              REAL(KIND=8) :: Q(1:32,-1:4,-1:4,-1:4,1:8)
              REAL(KIND=8) :: DQ(1:32,-1:4,-1:4,-1:4,1:8,1:3)
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE USLOPE
          END INTERFACE 
        END MODULE USLOPE__genmod
