        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:05 2020
        MODULE TSC_AMR__genmod
          INTERFACE 
            SUBROUTINE TSC_AMR(IND_CELL,IND_PART,IND_GRID_PART,X0,NG,NP,&
     &ILEVEL)
              INTEGER(KIND=4) :: IND_CELL(1:32)
              INTEGER(KIND=4) :: IND_PART(1:32)
              INTEGER(KIND=4) :: IND_GRID_PART(1:32)
              REAL(KIND=8) :: X0(1:32,1:3)
              INTEGER(KIND=4) :: NG
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: ILEVEL
            END SUBROUTINE TSC_AMR
          END INTERFACE 
        END MODULE TSC_AMR__genmod
