      SUBROUTINE JACOB (FN,F,X,XD,V,IO,JO,ABC,NR,NC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer NN,MM
      PARAMETER (NN=20, MM=10)      
      !DIMENSION X(*),XD(*),V(*),IO(*),JO(*),ABC(*)
      real(8) X(NN),XD(NN),V(NN),ABC(NN*NN)
      integer IO(NN),JO(NN)
      EXTERNAL FN,F
      LOGICAL FLAG, DIAGS
      CHARACTER::ANS
      REAL(8) FN,TDV
      DATA DEL,DMIN,TOLMIN,OKTOL /0.01d0, 0.05d0, 1D-6, 1D-4/
      
      !DIAGS= .TRUE.
      !PRINT '(1X,A,$)', 'DIAGNOSTICS ? (Y/N, "/"=N) '
      !READ(*,'(A)') ANS
      !IF (ANS== 'N' .or. ANS== 'n') DIAGS=.FALSE.
      DIAGS= .false.
      IJ= 1
      DO 40 J=1,NC
        DV= DMAX1( ABS( DEL*V(JO(J)) ), DMIN )
        DO 40 I=1,NR
          FLAG= .FALSE.
1         TOL= 0.1
          OLTOL= TOL
          TDV= DV 
          A2= 0.0
          A1= 0.0
          A0= 0.0
          B1= 0.0
          B0= 0.0
          D1= 0.0
          D0= 0.0
          IF (DIAGS .OR. FLAG) then
            WRITE(*,'(/1X,A8,I2,A1,I2,11X,A12,8X,A5)')
     &     'Element ',I,',',J, 'perturbation','slope'
          endif
          
          DO 20 K= 1,28 ! iterations on TDV
            A2= A1
            A1= A0
            B1= B0
            D1= D0
            A0= FN(F,XD,X,IO(I),JO(J),TDV)
            B0= DMIN1( ABS(A0), ABS(A1) )
            D0= ABS ( A0 - A1 )
            IF (DIAGS .OR. FLAG)then
              WRITE(*,'(20X,1P2E17.6)') TDV,A0
              goto 30
            endif
            
            IF(K .LE. 2) GO TO 20
            IF (A0 .EQ. A1 .AND. A1 .EQ. A2) THEN
              ANS2= A1
              GO TO 30
            END IF
            IF (A0 .EQ. 0.0) GO TO 25
10          IF( D0 .LE. TOL*B0 .AND. D1 .LE. TOL*B1) THEN
              ANS2= A1
              OLTOL= TOL
              IF(DIAGS .OR. FLAG)then
                WRITE(*,'(1X,A9,F8.7)') 'MET TOL= ',TOL
              endif
              IF (TOL .LE. TOLMIN) THEN
                GO TO 30
              ELSE
                TOL= 0.2*TOL
                GO TO 10
              END IF
            END IF
20          TDV= 0.6D0*TDV
      
25        IF (OLTOL .LE. OKTOL) THEN
            GO TO 30
          ELSE IF (.NOT. FLAG) THEN
            WRITE(*,'(/1X,A)')'NO CONVERGENCE *****'
            FLAG= .TRUE.
            GO TO 1
          ELSE
21          WRITE(*,'(1X,A,$)') 'Enter estimate : '
            READ(*,*,ERR=21) ANS2
            FLAG= .FALSE.
            GO TO 30
          END IF
30        ABC(IJ)= ANS2
          IF (DIAGS) THEN
            PRINT '(27X,A5,1PE13.6)','Ans= ',ANS2
            PAUSE 'Press "enter"'
          END IF
40        IJ= IJ+1
      RETURN 
      END

      !使用中心差分计算Df(I)/DX(J)
      FUNCTION FDX(F,XD,X,I,J,DDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer NN,MM
      PARAMETER (NN=20, MM=10)      
      REAL(8):: XD(NN), X(NN), FDX
      real(8):: T, DDX, XD1, XD2
      EXTERNAL F
      T =  X(J) 
      X(J)=  T - DDX 
      CALL F(X,XD)
      XD1 =  XD(I) 
      X(J)=  T + DDX 
      CALL F(X,XD)
      XD2 =  XD(I) 
      FDX = (XD2-XD1)/(DDX+DDX)
      X(J)=  T 
      RETURN
      END
      
      !使用中心差分计算Df(I)/DU(J)
      FUNCTION FDU(F,XD,X,I,J,DDU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer NN,MM
      PARAMETER (NN=20, MM=10)            
      PARAMETER (NIN=10)
      REAL(8)::XD(NN), X(NN), FDU
      COMMON/CONTROLS/U(NIN)
      real(8) T, DDU, XD1, XD2
      EXTERNAL F
      T =  U(J) 
      U(J)=  T - DDU 
      CALL F(X,XD)
      XD1 =  XD(I) 
      U(J)=  T + DDU 
      CALL F(X,XD)
      XD2 =  XD(I) 
      FDU = (XD2-XD1)/(DDU+DDU)
      U(J)=  T 
      RETURN
      END
      
      !使用中心差分计算DY(I)/DX(J)
      FUNCTION YDX(F,XD,X,I,J,DDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NOP=20)
      REAL(8) XD(*), X(*), YDX
      COMMON/OUTPUT/Y(NOP)
      REAL(8)::T, DDX, YD1, YD2
      EXTERNAL F
      T =  X(J) 
      X(J)=  T - DDX 
      CALL F(X,XD)
      YD1 =  Y(I) 
      X(J)=  T + DDX 
      CALL F(X,XD)
      YD2 =  Y(I) 
      YDX = (YD2-YD1)/(DDX+DDX)
      X(J)= T
      RETURN
      END
      
      !使用中心差分计算DY(I)/DU(J)
      FUNCTION YDU(F,XD,X,I,J,DDU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NIN=10, NOP=20)
      REAL(8)::XD(*), X(*), YDU
      COMMON/CONTROLS/U(NIN)
      COMMON/OUTPUT/Y(NOP)
      REAL(8)::T, DDU, YD1, YD2
      EXTERNAL F
      T =  U(J) 
      U(J)=  T - DDU 
      CALL F(X,XD)
      YD1 =  Y(I) 
      U(J)=  T + DDU 
      CALL F(X,XD)
      YD2 =  Y(I) 
      YDU = (YD2-YD1)/(DDU+DDU)
      U(J)= T
      RETURN
      END
      