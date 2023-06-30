        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:01 2020
        MODULE STELLAR_MOMENTUM__genmod
          INTERFACE 
            SUBROUTINE STELLAR_MOMENTUM(PIN,FLUX,DX,DY,DZ,DT,NGRID)
              REAL(KIND=8) :: PIN(1:32,-1:4,-1:4,-1:4)
              REAL(KIND=8) :: FLUX(1:32,1:3,1:3,1:3,1:8,1:3)
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DY
              REAL(KIND=8) :: DZ
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE STELLAR_MOMENTUM
          END INTERFACE 
        END MODULE STELLAR_MOMENTUM__genmod
