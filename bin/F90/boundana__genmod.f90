        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:03 2020
        MODULE BOUNDANA__genmod
          INTERFACE 
            SUBROUTINE BOUNDANA(X,U,IND_CELL,X_REF,IND_CELL_REF,DX,     &
     &IBOUND,NCELL)
              REAL(KIND=8) :: X(1:32,1:3)
              REAL(KIND=8) :: U(1:32,1:8)
              INTEGER(KIND=4) :: IND_CELL(1:32)
              REAL(KIND=8) :: X_REF(1:32,1:3)
              INTEGER(KIND=4) :: IND_CELL_REF(1:32)
              REAL(KIND=8) :: DX
              INTEGER(KIND=4) :: IBOUND
              INTEGER(KIND=4) :: NCELL
            END SUBROUTINE BOUNDANA
          END INTERFACE 
        END MODULE BOUNDANA__genmod
