      !给定状态变量X，根据6DOF方程计算X对时间的导数XD
!                       1   2     3    4   5   6  7 8 9 10       11       12     
!    X = STATE VECTOR - VT,ALPHA,BETA,PHI,THE,PSI,P,Q,R,XN(north),XE(east),H.
      SUBROUTINE F(X,XD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(*), XD(*), D(9), MASS
      logical::Ltrim
      COMMON/logicalvar/Ltrim
      COMMON/PARAM/XCG
      COMMON/CONTROLS/THTL,EL,AIL,RDR
      COMMON/OUTPUT/AN,ALAT,AX,QBAR,AMACH,Q,ALPHA,T
      PARAMETER (AXX=12874.847d0, AYY=75673.624d0,AZZ=85552.113d0,
     &           AXZ=1331.413d0)
      PARAMETER (AXZS=AXZ**2, XPQ=AXZ*(AXX-AYY+AZZ),
     &           GAM=AXX*AZZ-AXZ**2)
      PARAMETER (XQR= AZZ*(AZZ-AYY)+AXZS, ZPQ=(AXX-AYY)*AXX+AXZS)
      PARAMETER ( YPR= AZZ - AXX )
      PARAMETER ( GD=9.8066d0,mass=9298.644d0)
      COMMON /PBLK/ S,B,CBAR,XCGR,HX,G,RTOD
      DATA S,B,CBAR,XCGR,HX,G
     & /27.871d0,9.144d0,3.45d0,0.35d0,216.93d0,9.8066d0/
      !S planform area, ft^2->m^2
      !B span, ft -> m
      !CBAR mean aero chord, ft->m
      !XCGR reference center of gravity as a fraction of cbar
      
      RTOD=180.d0/acos(-1.d0)
C
C Assign state & control variables
C
      VT= X(1); ALPHA= X(2)*RTOD; BETA= X(3)*RTOD
      PHI=X(4); THETA= X(5); PSI= X(6)
      P= X(7); Q= X(8); R= X(9); ALT= X(12)
      U=X(1)*COS(X(2))*COS(X(3))
      V=X(1)          *SIN(X(3))
      W=X(1)*SIN(X(2))*COS(X(3))

      CBTA = COS(X(3));
C
C Air data computer and engine model
C
      CALL ADC(VT,ALT,AMACH,QBAR)
      !Ltrim=.true.时，时域响应时油门开度不变
      !Ltrim=.false.时，时域响应是推力不变
      if(Ltrim)then
        T= ENGINE_THRUST(THTL,ALT,AMACH)
      else
        T=T
      endif
      
! Look-up tables and component buildup

      call F16_AERO(CXT,CYT,CZT,CLT,CMT,CNT,
     $                    VT,ALPHA,BETA,P,Q,R,
     $                    EL,AIL,RDR,XCG)
      
      
!C Get ready for state equations


      STH = SIN(THETA); CTH= COS(THETA); SPH= SIN(PHI)
      CPH = COS(PHI) ; SPSI= SIN(PSI); CPSI= COS(PSI)
      QS = QBAR * S ; QSB= QS * B; RMQS= QS/MASS
      GCTH = GD * CTH ; QSPH= Q * SPH
      AY = RMQS*CYT ; AZ= RMQS * CZT

C Force equations

      UDOT = R*V - Q*W - GD*STH + (QS * CXT + T)/MASS
      VDOT = P*W - R*U + GCTH * SPH + AY
      WDOT = Q*U - P*V + GCTH * CPH + AZ
      DUM = (U*U + W*W)
      xd(1) = (U*UDOT + V*VDOT + W*WDOT)/VT
      xd(2) = (U*WDOT - W*UDOT) / DUM
      xd(3) = (VT*VDOT- V*XD(1)) * CBTA / DUM
      xd(13)= UDOT
      xd(14)= VDOT
      xd(15)= WDOT

C Kinematics

      xd(4) = P + (STH/CTH)*(QSPH + R*CPH)
      xd(5) = Q*CPH - R*SPH
      xd(6) = (QSPH + R*CPH)/CTH

C Moments

      ROLL = QSB*CLT
      PITCH = QS *CBAR*CMT
      YAW = QSB*CNT
      PQ = p*Q
      QR = Q*R
      QHX = Q*HX
      xd(7) = ( XPQ*PQ - XQR*QR + AZZ*ROLL + AXZ*(YAW + QHX) )/GAM
      xd(8) = ( YPR*P*R - AXZ*(P**2 - R**2) + PITCH - R*HX )/AYY
      xd(9) = ( ZPQ*PQ - XPQ*QR + AXZ*ROLL + AXX*(YAW + QHX) )/GAM

C Navigation

      T1= SPH * CPSI; T2= CPH * STH; T3= SPH * SPSI
      S1= CTH * CPSI; S2= CTH * SPSI; S3= T1 * STH - CPH * SPSI
      S4= T3 * STH + CPH * CPSI; S5= SPH * CTH; S6= T2*CPSI + T3
      S7= T2 * SPSI - T1; S8= CPH * CTH
C
      xd(10) = U * S1 + V * S3 + W * S6 ! North speed
      xd(11) = U * S2 + V * S4 + W * S7 ! East speed
      xd(12) = U * STH -V * S5 - W * S8 ! Vertical speed

C Outputs

      AN = -AZ/GD; ALAT= AY/GD
      RETURN
      END

      SUBROUTINE linearF(XD,AAA,BBB)
      use var
      IMPLICIT real(8) (A-H,O-Z)
      integer NX,NU!状态变量和控制变量的个数
      PARAMETER (NX=4, NU=2)      
      real(8) XD(NX)
      real(8) AAA(NX,NX),BBB(NX,NU)

      RTOD=180.d0/acos(-1.d0)
      do i=1,NX
        XD(i)=0.d0
        do j=1,NX
          XD(i)=XD(i)+AAA(i,j)*DeltaX(j)
        enddo
        do j=1,NU
          XD(i)=XD(i)+BBB(i,j)*DeltaU(j)
        enddo
      enddo
      RETURN
      END
    