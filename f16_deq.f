      SUBROUTINE F16_DEQ(XD,U,X,C,NI,NS,NC)
C
C  USAGE: CALL F16_DEQ(XD,U,X,C,NI,NS,NC)
C
C  DESCRIPTION:
C
C    COMPUTES STATE DERIVATIVES FOR THE F-16 MODEL.
C
C  INPUT:
C
C    U = INPUT VECTOR - THTL,EL,AIL,RDR.+ VXTURB,VYTERB,VZTERB
C    UPARM = ADD INPUT VECTOR - DELXCG, DELCX,DELCY,DELCZ,DELCL,DELCM,DELCN
C    X = STATE VECTOR - VT,ALPHA,BETA,PHI,THE,PSI,P,Q,R,XN(orth),XE(east),H.
C    C = VECTOR OF CONSTANTS - C(1) THROUGH C(9) = INERTIA TENSOR ELEMENTS.
C                              C(10) = AIRCRAFT MASS, SLUGS.
C                              C(11) = XCG, LONGITUDINAL C.G. LOCATION, 
C                                      DISTANCE NORMALIZED BY THE M.A.C.
C    NI = NUMBER OF INPUTS.
C    NS = NUMBER OF STATES.
C    NC = NUMBER OF ELEMENTS FOR INPUT VECTOR C.
C
C  OUTPUT
C
C    XD = STATE VECTOR TIME DERIVATIVE.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION XD(NS)
      DIMENSION U(NI),X(NS),C(NC)
C
      COMMON /PBLK/ S,B,CBAR,XCGR,HE,G,RTOD
C
      DATA S,B,CBAR,XCGR,HE,G,RTOD/
     $     300.,30.,11.32,0.35,160.,32.174,57.29578/
C
      VT1=X(1)
      ALPHA1=X(2)*RTOD
      BETA1=X(3)*RTOD
      PHI=X(4)
      THE=X(5)
      PSI=X(6)
      P=X(7)
      Q=X(8)
      R=X(9)
      ALT=X(12)
      POW=X(13)
      THTL=U(1)
      EL=U(2)
      AIL=U(3)
      RDR=U(4)
      RM=1.0/C(10)
      XCG=C(11)
C
C  TURBULENCE 
C
      VXTURB=U(5)
      VYTURB=U(6)
      VZTURB=U(7)
C
C  INCORPORATING GUST INTO AERO COEFFICIENTS TABLE LOOK UP
C
      CBTA=COS(X(3))
      VX1=VT1*COS(X(2))*CBTA
      VY1=VT1*SIN(X(3))
      VZ1=VT1*SIN(X(2))*CBTA
C      
C      
      VX=VX1+VXTURB
      VY=VY1+VYTURB
      VZ=VZ1+VZTURB      
C
      VT=SQRT(VX*VX+VY*VY+VZ*VZ)
      ALPHA=ATAN(VZ/VX)*RTOD
      BETA=ASIN(VY/VT)*RTOD
C
      CALL ATM(VT,ALT,RMACH,QBAR)
      CALL F16_ENGINE(THRUST,DPOW,POW,ALT,RMACH,THTL)
      XD(13)=DPOW
      CALL F16_AERO(CX,CY,CZ,CL,CM,CN,
     $              VT,ALPHA,BETA,P,Q,R,
     $              EL,AIL,RDR,XCG)

      WRITE(6,101) CX
 101  FORMAT(/,' CX = ',E13.6)
      WRITE(6,102) CY
 102  FORMAT(/,' CY = ',E13.6)
      WRITE(6,103) CZ
 103  FORMAT(/,' CZ = ',E13.6)
      WRITE(6,104) CL
 104  FORMAT(/,' CL = ',E13.6)
      WRITE(6,105) CM
 105  FORMAT(/,' CM = ',E13.6)
      WRITE(6,106) CN
 106  FORMAT(/,' CN = ',E13.6)
      CX=CX;
      CY=CY;
      CZ=CZ;
      CL=CL;
      CM=CM;
      CN=CN;
      STH=SIN(THE)
      CTH=COS(THE)
      SPH=SIN(PHI)
      CPH=COS(PHI)
      SPSI=SIN(PSI)
      CPSI=COS(PSI)
      QS=QBAR*S
      AY=RM*QS*CY
      AZ=RM*QS*CZ
C
C  FORCE EQUATIONS.
C
      VXDOT=R*VY - Q*VZ - G*STH + RM*(QS*CX + THRUST)
      VYDOT=P*VZ - R*VX + G*CTH*SPH + AY
      VZDOT=Q*VX - P*VY + G*CTH*CPH + AZ
      DEN=(VX*VX + VZ*VZ)
      XD(1)=(VX*VXDOT + VY*VYDOT + VZ*VZDOT)/VT
      XD(2)=(VX*VZDOT -VZ*VXDOT)/DEN
      XD(3)=(VT*VYDOT - VY*XD(1))*CBTA/DEN
C
C  KINEMATICS.
C     
      XD(5)=Q*CPH - R*SPH
      XD(6)=(Q*SPH + R*CPH)/CTH
	  XD(4)=P + STH*XD(6)
C
C  MOMENT EQUATIONS.
C
      XD(7)=(C(1)*R + C(2)*P + C(4)*HE)*Q + QS*B*(C(3)*CL + C(4)*CN)
      XD(8)=(C(5)*P - C(7)*HE)*R + C(6)*(R*R-P*P) + C(7)*QS*CBAR*CM
      XD(9)=(C(8)*P - C(2)*R + C(9)*HE)*Q + QS*B*(C(4)*CL + C(9)*CN)
C
C  NAVIGATION EQUATIONS.
C
      XD(10)=VX*CTH*CPSI + VY*(SPH*STH*CPSI-CPH*SPSI)
     $       + VZ*(SPH*SPSI+CPH*STH*CPSI)
      XD(11)=VX*CTH*SPSI + VY*(CPH*CPSI+SPH*STH*SPSI)
     $       + VZ*(CPH*STH*SPSI-SPH*CPSI)
      XD(12)=VX*STH - VY*SPH*CTH - VZ*CPH*CTH
C
      RETURN
      END
C
C
      SUBROUTINE ATM(VT,ALT,RMACH,QBAR)
C
C  This subroutine computes properties of the standard atmosphere.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      RO = 0.002377
      TFAC = 1.0 - ALT*0.703E-05
      T = 519.0*TFAC
      IF(ALT.GE.35000.0) T=390.0
      RHO = RO*(TFAC**4.14)
      RMACH = VT/SQRT(1.4*1716.3*T)
      QBAR = 0.5*RHO*VT*VT
C      PS = 1715.0*RHO*T
C
      RETURN
      END
C
C
      SUBROUTINE F16_AERO(CX,CY,CZ,CL,CM,CN,
     $                    VT,ALPHA,BETA,P,Q,R,
     $                    EL,AIL,RDR,XCG)
C
C  This subroutine computes the body axis nondimensional aerodynamic
C  coefficients for the F-16 based on wind tunnel data
C  from NASA TP 1538, December 1979.  
C
C  INPUTS:
C
C    VT = AIRSPEED, FT/SEC.
C    ALPHA = ANGLE OF ATTACK, DEG.  ( -10 <= ALPHA <= 45 )
C    BETA = SIDESLIP ANGLE, DEG.  ( -30 <= BETA <= 30 )
C    P,Q,R = BODY AXIS ANGULAR VELOCITIES IN ROLL, PITCH, AND YAW, RAD/SEC.
C    EL = ELEVATOR DEFLECTION, DEG.  ( -25 <= EL <= 25 )
C    AIL = AILERON DEFLECTION, DEG.  ( -21.5 <= AIL <= 21.5 )
C    RDR = RUDDER DEFLECTION, DEG.  ( -30 <= RDR <= 30 )
C    XCG = LONGITUDINAL C.G. LOCATION, DISTANCE NORMALIZED BY THE M.A.C.
C
C
C  OUTPUTS:
C
C    CX,CY,CZ = BODY AXIS NONDIMENSIONAL AERODYNAMIC FORCE COEFFICIENTS.
C    CL,CM,CN = BODY AXIS NONDIMENSIONAL AERODYNAMIC MOMENT COEFFICIENTS.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION D(9)
C
      COMMON /PBLK/ S,B,CBAR,XCGR,HE,G,RTOD
C
C
C  BASIC FLOW ANGLE AND CONTROL SURFACE AERODYNAMICS.
C
      CX=CX_AERO(ALPHA,EL)
      CY=CY_AERO(BETA,AIL,RDR)
      CZ=CZ_AERO(ALPHA,BETA,EL)
      DAIL=AIL/20.0
      DRDR=RDR/30.0
      CL=CL_AERO(ALPHA,BETA)
      DCLDA=DLDA(ALPHA,BETA)
      DCLDR=DLDR(ALPHA,BETA)
      CL=CL + DCLDA*DAIL + DCLDR*DRDR
      CM=CM_AERO(ALPHA,EL)
      CN=CN_AERO(ALPHA,BETA)
      DCNDA=DNDA(ALPHA,BETA)
      DCNDR=DNDR(ALPHA,BETA)
      CN=CN + DCNDA*DAIL + DCNDR*DRDR
C
C  ADD DAMPING TERMS.
C
      CALL DAMP(ALPHA,D)
      CQ=0.5*Q*CBAR/VT
      B2V=0.5*B/VT
      CX=CX + CQ*D(1)
      CY=CY + B2V*(D(2)*R + D(3)*P)
      CZ=CZ + CQ*D(4)
      CL=CL + B2V*(D(5)*R + D(6)*P)
      CM=CM + CQ*D(7) + CZ*(XCGR-XCG)
      CN=CN + B2V*(D(8)*R + D(9)*P) - CY*(XCGR-XCG)*CBAR/B
C
      RETURN
      END
C
C
      SUBROUTINE DAMP(ALPHA, D)
C
C  This subroutine computes the damping coefficients
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C
C  D = [CXq, CYr, CYp, CZq, Clr, Clp, Cmq, Cnr, Cnp]
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(-2:9,9),D(9)
C
      DATA A/  -.267,-.110,.308,1.34,2.08,2.91,2.76,2.05,
     $         1.50,1.49,1.83,1.21,
     $          .882,.852,.876,.958,.962,.974,.819,.483,
     $          .590,1.21,-.493,-1.04,
     $         -.108,-.108,-.188,.110,.258,.226,.344,.362,
     $          .611,.529,.298,-.227,
     $        -8.80,-25.8,-28.9,-31.4,-31.2,-30.7,-27.7,
     $       -28.2,-29.0,-29.8,-38.3,-35.3,
     $         -.126,-.026,.063,.113,.208,.230,.319,.437,
     $          .680,.100,.447,-.330,
     $         -.360,-.359,-.443,-.420,-.383,-.375,-.329,-.294,
     $         -.230,-.210,-.120,-.100,
     $        -7.21,-5.40,-5.23,-5.26,-6.11,-6.64,-5.69,
     $        -6.00,-6.20,-6.40,-6.60,-6.00,
     $         -.380,-.363,-.378,-.386,-.370,-.453,-.550,
     $         -.582,-.595,-.637,-1.02,-.840,
     $          .061,.052,.052,-.012,-.013,-.024,.050,.150,
     $          .130,.158,.240,.150 /
C
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      DO 10  I=1,9
        D(I) = A(K,I) + ABS(DA)*(A(L,I)-A(K,I))
 10   CONTINUE
C
      RETURN
      END
C
C
      FUNCTION CX_AERO(ALPHA,EL)
C
C  This function computes the X body axis aerodynamic force coefficient
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  EL = ELEVATOR DEFLECTION, DEG  ( -25 <= EL <= 25 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(-2:9,-2:2)
C
      DATA A/  -.099,-.081,-.081,-.063,-.025,.044,
     $          .097,.113,.145,.167,.174,.166,
     $         -.048,-.038,-.040,-.021,.016,.083,
     $          .127,.137,.162,.177,.179,.167,
     $         -.022,-.020,-.021,-.004,.032,.094,
     $          .128,.130,.154,.161,.155,.138,
     $         -.040,-.038,-.039,-.025,.006,.062,
     $          .087,.085,.100,.110,.104,.091,
     $         -.083,-.073,-.076,-.072,-.046,.012,
     $          .024,.025,.043,.053,.047,.040 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = EL/12.0
      M = INT(S)
      IF(M.LE.-2) M=-1
      IF(M.GE.2) M=1
      DE = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DE))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      CX_AERO = V + (W-V)*ABS(DE)
C
      RETURN
      END
C
C
      FUNCTION CY_AERO(BETA,AIL,RDR)
C
C  This function computes the Y body axis aerodynamic force coefficient
C  for the F-16 aerodynamic model.  
C
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C  AIL = AILERON DEFLECTION, DEG  ( -21.5 <= AIL <= 21.5 )
C  RDR = RUDDER DEFLECTION, DEG  ( -30 <= RDR <= 30 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      CY_AERO = -.02*BETA + .021*(AIL/20.0) + .086*(RDR/30.0)
C
      RETURN
      END
C
C
      FUNCTION CZ_AERO(ALPHA,BETA,EL)
C
C  This function computes the Z body axis aerodynamic force coefficient
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C  EL = ELEVATOR DEFLECTION, DEG  ( -25 <= EL <= 25 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(-2:9)
C
      DATA A /  .770,.241,-.100,-.416,-.731,-1.053,
     $        -1.366,-1.646,-1.917,-2.120,-2.248,-2.229 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = A(K) + ABS(DA)*(A(L)-A(K))
      CZ_AERO = S*(1-(BETA/57.3)**2) - .19*(EL/25.0)
C
      RETURN
      END
C
C
      FUNCTION CM_AERO(ALPHA,EL)
C
C  This function computes the Y body axis aerodynamic moment coefficient
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  EL = ELEVATOR DEFLECTION, DEG  ( -25 <= EL <= 25 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(-2:9,-2:2)
C
      DATA A/  .205,.168,.186,.196,.213,.251,
     $         .245,.238,.252,.231,.198,.192,
     $         .081,.077,.107,.110,.110,.141,
     $         .127,.119,.133,.108,.081,.093,
     $        -.046,-.020,-.009,-.005,-.006,.010,
     $         .006,-.001,.014,.000,-.013,.032,
     $        -.174,-.145,-.121,-.127,-.129,-.102,
     $        -.097,-.113,-.087,-.084,-.069,-.006,
     $        -.259,-.202,-.184,-.193,-.199,-.150,
     $        -.160,-.167,-.104,-.076,-.041,-.005/
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = EL/12.0
      M = INT(S)
      IF(M.LE.-2) M=-1
      IF(M.GE.2) M=1
      DE = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DE))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      CM_AERO = V + (W-V)*ABS(DE)
C
      RETURN
      END
C
C
      FUNCTION CL_AERO(ALPHA,BETA)
C
C  This function computes the X body axis aerodynamic moment coefficient
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(-2:9,0:6)
C
      DATA A/  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     $         -.001,-.004,-.008,-.012,-.016,-.022,
     $         -.022,-.021,-.015,-.008,-.013,-.015,
     $         -.003,-.009,-.017,-.024,-.030,-.041,
     $         -.045,-.040,-.016,-.002,-.010,-.019,
     $         -.001,-.010,-.020,-.030,-.039,-.054,
     $         -.057,-.054,-.023,-.006,-.014,-.027,
     $          .000,-.010,-.022,-.034,-.047,-.060,
     $         -.069,-.067,-.033,-.036,-.035,-.035,
     $          .007,-.010,-.023,-.034,-.049,-.063,
     $         -.081,-.079,-.060,-.058,-.062,-.059,
     $          .009,-.011,-.023,-.037,-.050,-.068,
     $         -.089,-.088,-.091,-.076,-.077,-.076 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = .2*ABS(BETA)
      M = INT(S)
      IF(M.EQ.0) M=1
      IF(M.GE.6) M=5
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      CL_AERO= (V + (W-V)*ABS(DB))*SIGN(1.0,BETA)
C
      RETURN
      END
C
C
      FUNCTION CN_AERO(ALPHA,BETA)
C
C  This function computes the Z body axis aerodynamic moment coefficient
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(-2:9,0:6)
C
      DATA A/  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     $          .018,.019,.018,.019,.019,.018,
     $          .013,.007,.004,-.014,-.017,-.033,
     $          .038,.042,.042,.042,.043,.039,
     $          .030,.017,.004,-.035,-.047,-.057,
     $          .056,.057,.059,.058,.058,.053,
     $          .032,.012,.002,-.046,-.071,-.073,
     $          .064,.077,.076,.074,.073,.057,
     $          .029,.007,.012,-.034,-.065,-.041,
     $          .074,.086,.093,.089,.080,.062,
     $          .049,.022,.028,-.012,-.002,-.013,
     $          .079,.090,.106,.106,.096,.080,
     $          .068,.030,.064,.015,.011,-.001 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = .2*ABS(BETA)
      M = INT(S)
      IF(M.EQ.0) M=1
      IF(M.GE.6) M=5
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      CN_AERO= (V + (W-V)*ABS(DB))*SIGN(1.0,BETA)
C
      RETURN
      END
C
C
      FUNCTION DLDA(ALPHA,BETA)
C
C  This function computes the rolling moment due to aileron deflection
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(-2:9,-3:3)
C
      DATA A/  -.041,-.052,-.053,-.056,-.050,-.056,
     $         -.082,-.059,-.042,-.038,-.027,-.017,
     $         -.041,-.053,-.053,-.053,-.050,-.051,
     $         -.066,-.043,-.038,-.027,-.023,-.016,
     $         -.042,-.053,-.052,-.051,-.049,-.049,
     $         -.043,-.035,-.026,-.016,-.018,-.014,
     $         -.040,-.052,-.051,-.052,-.048,-.048,
     $         -.042,-.037,-.031,-.026,-.017,-.012,
     $         -.043,-.049,-.048,-.049,-.043,-.042,
     $         -.042,-.036,-.025,-.021,-.016,-.011,
     $         -.044,-.048,-.048,-.047,-.042,-.041,
     $         -.020,-.028,-.013,-.014,-.011,-.010,
     $         -.043,-.049,-.047,-.045,-.042,-.037,
     $         -.003,-.013,-.010,-.003,-.007,-.008 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = 0.1*BETA
      M = INT(S)
      IF(M.EQ.-3) M=-2
      IF(M.GE.3) M=2
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      DLDA = V + (W-V)*ABS(DB)
C
      RETURN
      END
C
C
      FUNCTION DLDR(ALPHA,BETA)
C
C  This function computes the rolling moment due to rudder deflection
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(-2:9,-3:3)
C
      DATA A/  .005,.017,.014,.010,-.005,.009,
     $         .019,.005,.000,-.005,-.011,.008,
     $         .007,.016,.014,.014,.013,.009,
     $         .012,.005,.000,.004,.009,.007,
     $         .013,.013,.011,.012,.011,.009,
     $         .008,.005,.000,.005,.003,.005,
     $         .018,.015,.015,.014,.014,.014,
     $         .014,.015,.013,.011,.006,.001,
     $         .015,.014,.013,.013,.012,.011,
     $         .011,.010,.008,.008,.007,.003,
     $         .021,.011,.010,.011,.010,.009,
     $         .008,.010,.006,.005,.000,.001,
     $         .023,.010,.011,.011,.011,.010,
     $         .008,.010,.006,.014,.020,.000 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = 0.1*BETA
      M = INT(S)
      IF(M.EQ.-3) M=-2
      IF(M.GE.3) M=2
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      DLDR = V + (W-V)*ABS(DB)
C
      RETURN
      END
C
C
      FUNCTION DNDA(ALPHA,BETA)
C
C  This function computes the yawing moment due to aileron deflection
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(-2:9,-3:3)
C
      DATA A/  .001,-.027,-.017,-.013,-.012,-.016,
     $         .001,.017,.011,.017,.008,.016,
     $         .002,-.014,-.016,-.016,-.014,-.019,
     $        -.021,.002,.012,.016,.015,.011,
     $        -.006,-.008,-.006,-.006,-.005,-.008,
     $        -.005,.007,.004,.007,.006,.006,
     $        -.011,-.011,-.010,-.009,-.008,-.006,
     $         .000,.004,.007,.010,.004,.010,
     $        -.015,-.015,-.014,-.012,-.011,-.008,
     $        -.002,.002,.006,.012,.011,.011,
     $        -.024,-.010,-.004,-.002,-.001,.003,
     $         .014,.006,-.001,.004,.004,.006,
     $        -.022,.002,-.003,-.005,-.003,-.001,
     $        -.009,-.009,-.001,.003,-.002,.001 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = 0.1*BETA
      M = INT(S)
      IF(M.EQ.-3) M=-2
      IF(M.GE.3) M=2
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      DNDA = V + (W-V)*ABS(DB)
C
      RETURN
      END
C
C
      FUNCTION DNDR(ALPHA,BETA)
C
C  This function computes the yawing moment due to rudder deflection
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(-2:9,-3:3)
C
      DATA A/  -.018,-.052,-.052,-.052,-.054,-.049,
     $         -.059,-.051,-.030,-.037,-.026,-.013,
     $         -.028,-.051,-.043,-.046,-.045,-.049,
     $         -.057,-.052,-.030,-.033,-.030,-.008,
     $         -.037,-.041,-.038,-.040,-.040,-.038,
     $         -.037,-.030,-.027,-.024,-.019,-.013,
     $         -.048,-.045,-.045,-.045,-.044,-.045,
     $         -.047,-.048,-.049,-.045,-.033,-.016,
     $         -.043,-.044,-.041,-.041,-.040,-.038,
     $         -.034,-.035,-.035,-.029,-.022,-.009,
     $         -.052,-.034,-.036,-.036,-.035,-.028,
     $         -.024,-.023,-.020,-.016,-.010,-.014,
     $         -.062,-.034,-.027,-.028,-.027,-.027,
     $         -.023,-.023,-.019,-.009,-.025,-.010 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = 0.1*BETA
      M = INT(S)
      IF(M.EQ.-3) M=-2
      IF(M.GE.3) M=2
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      DNDR = V + (W-V)*ABS(DB)
C
      RETURN
      END
C
C
      SUBROUTINE F16_ENGINE(THRUST,DPOW,POW,ALT,RMACH,THTL)
C
C  This subroutine computes the engine thrust and the time derivative
C  of the engine power state for the F-16.
C
C  INPUTS:
C
C    POW = ENGINE POWER LEVEL, PERCENT.  ( 0 <= POW <= 100. )
C    ALT = ALTITUDE, FT.        ( 0 <= ALT <= 50000. )
C    RMACH = MACH NUMBER.       ( 0 <= RMACH <= 1.0 )
C    THTL = THROTTLE SETTING.   ( 0 <= THTL <= 1.0 )
C
C
C  OUTPUTS:
C
C    THRUST = ENGINE THRUST, LBF.  
C    DPOW = TIME DERIVATIVE OF THE ENGINE POWER LEVEL, PERCENT/SEC.
C
C
C  COMPUTE ENGINE THRUST.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      THRUST=ENGINE_THRUST(POW,ALT,RMACH)
C
C  COMPUTE COMMANDED POWER LEVEL AND POWER LEVEL TIME DERIVATIVE.
C
      CPOW=TGEAR(THTL)
      DPOW=PDOT(POW,CPOW)
C
      RETURN
      END
C
C
      FUNCTION ENGINE_THRUST(POW,ALT,RMACH)
C
C  This function computes the thrust 
C  for the F-16 model.
C
C  POW = ENGINE POWER LEVEL, PERCENT.  ( 0 <= POW <= 100. )
C  ALT = ALTITUDE, FT.        ( 0 <= ALT <= 50000. )
C  RMACH = MACH NUMBER.       ( 0 <= RMACH <= 1.0 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION A(0:5,0:5), B(0:5,0:5), C(0:5,0:5)
C
C  IDLE POWER DATA.
C
      DATA A/  1060.,  670.,  880., 1140., 1500., 1860.,
     $          635.,  425.,  690., 1010., 1330., 1700.,
     $           60.,   25.,  345.,  755., 1130., 1525.,
     $        -1020., -710., -300.,  350.,  910., 1360.,
     $        -2700.,-1900.,-1300., -247.,  600., 1100.,
     $        -3600.,-1400., -595., -342., -200.,  700./
C
C  MIL POWER DATA.
C
      DATA B/ 12680., 9150., 6200., 3950., 2450., 1400.,
     $        12680., 9150., 6313., 4040., 2470., 1400.,
     $        12610., 9312., 6610., 4290., 2600., 1560.,
     $        12640., 9839., 7090., 4660., 2840., 1660.,
     $        12390.,10176., 7750., 5320., 3250., 1930.,
     $        11680., 9848., 8050., 6100., 3800., 2310./
C
C  MAX POWER DATA.
C
      DATA C/ 20000.,15000.,10800., 7000., 4000., 2500.,
     $        21420.,15700.,11225., 7323., 4435., 2600.,
     $        22700.,16860.,12250., 8154., 5000., 2835.,
     $        24240.,18910.,13760., 9285., 5700., 3215.,
     $        26070.,21075.,15975.,11115., 6860., 3950.,
     $        28886.,23319.,18300.,13484., 8642., 5057./
C
C
C  ROW INDEX FOR ALTITUDE.
C
      H=0.0001*ALT
      I=INT(H)
      IF (I.GE.5) I=4
      DH=H-FLOAT(I)
C
C  COLUMN INDEX FOR MACH NUMBER.
C
      RM=5.*RMACH
      M=INT(RM)
      IF (M.GE.5) M=4
      DM=RM-FLOAT(M)
      CDH=1.0-DH
C
C  COMPUTE MIL THRUST.
C
C  ALTITUDE INTERPOLATION.
      S= B(I,M)*CDH + B(I+1,M)*DH
      T= B(I,M+1)*CDH + B(I+1,M+1)*DH
C  MACH NUMBER INTERPOLATION.
      TMIL= S + (T-S)*DM
C
C  INTERPOLATE WITH IDLE OR MAX THRUST, DEPENDING ON POWER LEVEL COMMAND.
C
      IF (POW.LT.50.0) THEN
C
C  COMPUTE IDLE THRUST.
C
C  ALTITUDE INTERPOLATION.
        S= A(I,M)*CDH + A(I+1,M)*DH
        T= A(I,M+1)*CDH + A(I+1,M+1)*DH
C  MACH NUMBER INTERPOLATION.
        TIDL= S + (T-S)*DM
        ENGINE_THRUST= TIDL + (TMIL-TIDL)*POW/50.0
      ELSE
C
C  COMPUTE MAX THRUST.
C
C  ALTITUDE INTERPOLATION.
        S= C(I,M)*CDH + C(I+1,M)*DH
        T= C(I,M+1)*CDH + C(I+1,M+1)*DH
C  MACH NUMBER INTERPOLATION.
        TMAX= S + (T-S)*DM
        ENGINE_THRUST= TMIL + (TMAX-TMIL)*(POW-50.0)*0.02
      END IF
C
      RETURN
      END
C
C
      FUNCTION TGEAR(THTL)
C
C  This function computes the engine
C  power level command, POW, for an input throttle setting, THTL,
C  for the F-16 engine model.
C
C  THTL = THROTTLE SETTING.  ( 0 <= THTL <= 1.0 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      IF (THTL.LE.0.77) THEN
        TGEAR = 64.94*THTL
      ELSE
        TGEAR = 217.38*THTL-117.38
      END IF
      RETURN
      END
C
C
      FUNCTION PDOT(POW,CPOW)
C
C  This function computes the rate of change in engine power level
C  using a first order lag as a function of actual power, POW, 
C  and commanded power, CPOW, for the F-16 engine model.
C
C  POW =  ENGINE POWER LEVEL, PERCENT.  ( 0 <= POW <= 100. )
C  CPOW = COMMANDED ENGINE POWER LEVEL, PERCENT.  ( 0 <= CPOW <= 100. )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      IF (CPOW.GE.50.0) THEN
        IF (POW.GE.50.0) THEN
          TPOW=CPOW
          T=5.0
        ELSE
          TPOW=60.0
          T=RTAU(TPOW-POW)
        END IF
      ELSE
        IF (POW.GE.50.0) THEN
          TPOW=40.0
          T=5.0
        ELSE
          TPOW=CPOW
          T=RTAU(TPOW-POW)
        END IF
      END IF
      PDOT=T*(TPOW-POW)
      RETURN
      END
C
C
      FUNCTION RTAU(DP)
C
C  This function computes the thrust lag reciprocal time constant
C  for the F-16 engine model.
C
C  DP = CHANGE IN POWER LEVEL, PERCENT  ( 0 <= DP <= 100. )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      IF (DP.LE.25.0) THEN
        RTAU=1.0
      ELSE IF (DP.GE.50.0) THEN
        RTAU = 0.1
      ELSE
        RTAU = 1.9 - 0.036*DP
      END IF
      RETURN
      END