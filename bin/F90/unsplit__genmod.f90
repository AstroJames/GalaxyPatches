        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:02 2020
        MODULE UNSPLIT__genmod
          INTERFACE 
            SUBROUTINE UNSPLIT(UIN,GRAVIN,PIN,FLUX,TMP,DX,DY,DZ,DT,NGRID&
     &)
              REAL(KIND=8) :: UIN(1:32,-1:4,-1:4,-1:4,1:8)
              REAL(KIND=8) :: GRAVIN(1:32,-1:4,-1:4,-1:4,1:3)
              REAL(KIND=8) :: PIN(1:32,-1:4,-1:4,-1:4)
              REAL(KIND=8) :: FLUX(1:32,1:3,1:3,1:3,1:8,1:3)
              REAL(KIND=8) :: TMP(1:32,1:3,1:3,1:3,1:2,1:3)
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DY
              REAL(KIND=8) :: DZ
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE UNSPLIT
          END INTERFACE 
        END MODULE UNSPLIT__genmod
