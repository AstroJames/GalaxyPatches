        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:01 2020
        MODULE TRACEX__genmod
          INTERFACE 
            SUBROUTINE TRACEX(Q,DQ,C,QM,QP,DX,DT,NGRID)
              REAL(KIND=8) :: Q(1:32,-1:4,-1:4,-1:4,1:8)
              REAL(KIND=8) :: DQ(1:32,-1:4,-1:4,-1:4,1:8,1:3)
              REAL(KIND=8) :: C(1:32,-1:4,-1:4,-1:4)
              REAL(KIND=8) :: QM(1:32,-1:4,-1:4,-1:4,1:8,1:3)
              REAL(KIND=8) :: QP(1:32,-1:4,-1:4,-1:4,1:8,1:3)
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE TRACEX
          END INTERFACE 
        END MODULE TRACEX__genmod
