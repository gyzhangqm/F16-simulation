      SUBROUTINE TRIMMER (NV, COST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=20, MM=10)
      real(8) m2ft
      EXTERNAL COST
      CHARACTER*1 ANS
      DIMENSION S(6), DS(6)
      COMMON/ STATE/ X(NN)
      COMMON/ DSTATE/ XD(NN)
      COMMON/ CONTROLS/ U(MM)
      COMMON/ OUTPUT/ AN,AY,AX,QBAR,AMACH ! common to aircraft
      
      RTOD=180.d0/acos(-1.d0)
      m2ft=3.28084
      
      S(1) = U(1)
      S(2) = U(2)
      S(3) = X(2)
      DS(1) = 0.2d0
      DS(2) = 1.d0
      DS(3) = 0.02d0
      if(nv>3) then
        S(4)= U(3)
        S(5)= U(4)
        S(6)= X(3)
        DS(4)= 1.0d0
        DS(5)= 1.0d0
        DS(6)= 0.02d0
      endif
      
      NC= 5000
!      WRITE(*,'(1X,A,$)')'Reqd. # of trim iterations (def. = 1000) : '
!      READ(*,*,ERR=20) NC
      SIGMA = -1.0d0
      CALL SMPLX(COST,NV,S,DS,SIGMA,NC,F0,FFIN)
      FFIN = COST(S)
!      IF (NV .GT. 3) THEN
!        WRITE(*,'(/11X,A)')' Throttle    Elevator    Ailerons    Rudder'
!        WRITE(*,'(9X,4(1PE10.2,3X),/)') U(1), U(2), U(3), U(4)
!        WRITE(*,99)'Angle of attack',RTOD*X(2),
!     &             'Sideslip angle',RTOD*X(3)
!        WRITE(*,99) 'Pitch angle', RTOD*X(5), 'Bank angle', RTOD*X(4)
!        WRITE(*,99) 'Normal acceleration', AN, 'Lateral acceln', AY
!        WRITE(*,99) 'Dynamic pressure', QBAR, 'Mach number', AMACH
!      ELSE
!        WRITE(*,'(/1X,A)')'  Throttle    Elevator    Alpha       Pitch(��)'
!        WRITE(*,'(1X,4(1PE10.2,3X))')U(1),U(2),X(2)*RTOD,X(3)*RTOD
!        WRITE(*,'(/1X,A)')'Normal acceleration Dynamic Pressure Mach '
!        WRITE(*,'(5X,3(1PE10.2,7X))') AN,QBAR,AMACH
!      END IF

      WRITE(*,99)'Initial cost function ',F0,'Final cost function',FFIN
99    FORMAT(2(1X,A22,1PE10.2))
!40    WRITE(*,'(/1X,A,$)') 'More Iterations ? (def= Y) : '
!      READ(*,'(A)',ERR= 40) ANS
!      IF (ANS .EQ. 'Y'.OR. ANS .EQ. 'y') GO TO 10
!      IF (ANS .EQ. 'N'.OR. ANS .EQ. 'n') RETURN
!      GO TO 40
      

      !write(*,*) X(1:12)
      !write(*,*) "XD:"
      !write(*,*) XD(1:12)

      END

      FUNCTION COST (S) ! F16 cost function (see text)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=20)
      DOUBLE PRECISION S(*), XD(NN),X(NN)
      COMMON/STATE/X 
      COMMON/DSTATE/XD 
      COMMON/CONTROLS/THTL,EL,AIL,RDR ! to aircraft
      THTL = S(1)
      EL = S(2)
      X(2)= S(3)
      AIL = S(4)
      RDR = S(5)
      X(3) = S(6)
      !X(13)= TGEAR(THTL)
      CALL CONSTR (X)
      CALL F (X,XD)
      COST = XD(1)**2 + 100.0*( XD(2)**2 + XD(3)**2 )
     &     + 10.0*( XD(7)**2 + XD(8)**2 + XD(9)**2 )
      RETURN
      END

      SUBROUTINE CONSTR (X) ! used by COST, to apply constraints
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(*)
      LOGICAL COORD, STAB
      COMMON/CONSTRNT/RADGAM,SINGAM,RR,PR,TR,PHI,CPHI,SPHI,COORD,STAB

      CALPH = COS(X(2))
      SALPH = SIN(X(2))
      CBETA = COS(X(3))
      SBETA = SIN(X(3))

        X(3)= 0.d0
        X(4)= PHI
        
        D = X(2)
        IF(PHI .NE. 0.0) D = -X(2) ! inverted
        IF( SINGAM .NE. 0.0 ) THEN ! climbing
          SGOCB = SINGAM / CBETA
          X(5)= D + ATAN( SGOCB/SQRT(1.0-SGOCB*SGOCB)) ! roc constraint
        ELSE
          X(5) = D 
        END IF
        X(7)= RR
        X(8)= PR
        X(9) = 0.d0 ! body-axis roll
      RETURN
      END

      SUBROUTINE SMPLX(FX,N,X,DX,SD,M,Y0,YL)
C This simplex algorithm minimizes FX(X), where X is (Nx1).
C DX contains the initial perturbations in X. SD should be set according
C to the tolerance required; when SD<0 the algorithm calls FX M times
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL(8) X(*), DX(*)
      DIMENSION XX(32), XC(32), Y(33), V(32,32)

      NV=N+1
      DO I=1,N
        do J=1,NV
          V(I,J)=X(I) ! i娑撳搫娼楅弽鍥风礉j娑撹櫣顑噅娑擃亞�?
        enddo
        V(I,I+1)=X(I)+DX(I)
      enddo

      Y0=FX(X)
      Y(1)=Y0
      DO 3 J=2,NV
    3   Y(J)=FX(V(1,J))
      K=NV
    4 YH=Y(1)
      YL=Y(1)
      NH=1
      NL=1
      DO 5 J=2,NV
        IF(Y(J).GT.YH) THEN
          YH=Y(J)
          NH=J
        ELSEIF(Y(J).LT.YL) THEN
          YL=Y(J)
          NL=J
        ENDIF
    5 CONTINUE
      YB=Y(1)
      DO 6 J=2,NV
    6 YB=YB+Y(J)
      YB=YB/NV
      D=0.0
      DO 7 J=1,NV
    7 D=D+(Y(J)-YB)**2
      SDA=SQRT(D/NV)
      IF((K.GE.M).OR.(SDA.LE.SD)) THEN
        SD=SDA
        M=K
        YL=Y(NL)
        DO 8 I=1,N
    8     X(I)=V(I,NL)
        RETURN
      END IF
      DO 10 I=1,N
        XC(I)=0.0
        DO 9 J=1,NV
    9     IF(J.NE.NH) XC(I)=XC(I)+V(I,J)
   10     XC(I)=XC(I) /N
        DO 11 I=1,N
   11     X(I)=XC(I)+XC(I)-V(I,NH)
        K=K+1
        YR=FX(X)
        IF(YR.LT.YL) THEN
          DO 12 I=1,N
   12       XX(I)=X(I)+X(I)-XC(I)
          K=K+1
          YE=FX(XX)
          IF(YE.LT.YR) THEN
            Y(NH)=YE
            DO 13 I=1,N
   13         V(I,NH)=XX(I)
          ELSE
            Y(NH)=YR
            DO 14 I=1,N
   14         V(I,NH)=X(I)
          END IF
          GOTO 4
        ENDIF
        Y2=Y(NL)
        DO 15 J=1,NV
   15     IF((J.NE.NL).AND.(J.NE.NH).AND.(Y(J).GT.Y2)) Y2=Y(J)
            IF(YR.LT.YH) THEN
              Y(NH)=YR
              DO 16 I=1,N
   16           V(I,NH)=X(I)
              IF(YR.LT.Y2) GO TO 4
            ENDIF
        DO 17 I=1,N
   17     XX(I)=0.5*(V(I,NH)+XC(I))
        K=K+1
        YC=FX(XX)
        IF(YC.LT.YH) THEN
          Y(NH)=YC
          DO 18 I=1,N
   18       V(I,NH)=XX(I)
        ELSE
          DO 20 J=1,NV
            DO 19 I=1,N
   19         V(I,J)=0.5*(V(I,J)+V(I,NL))
   20      IF(J.NE.NL) Y(J)=FX(V(1,J))
             K=K+N
           ENDIF
        GO TO 4
      END