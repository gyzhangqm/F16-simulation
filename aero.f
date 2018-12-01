      SUBROUTINE F16_AERO(CX,CY,CZ,CL,CM,CN,
     $                    VT,ALPHA,BETA,P,Q,R,
     $                    EL,AIL,RDR,XCG)
C
C  This subroutine computes the body axis nondimensional aerodynamic
C  coefficients for the F-16 based on wind tunnel data
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
!   dail  = ail/20.0;   aileron normalized against max angle 
!   The aileron was normalized using 20.0 but the NASA report and
!   S&L both have 21.5 deg. as maximum deflection. 
!   rudder normalized against max angle 
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
     $1.50,1.49,1.83,1.21,
     $.882,.852,.876,.958,.962,.974,.819,.483,
     $.590,1.21,-.493,-1.04,
     $-.108,-.108,-.188,.110,.258,.226,.344,.362,
     $.611,.529,.298,-.227,
     $-8.80,-25.8,-28.9,-31.4,-31.2,-30.7,-27.7,
     $-28.2,-29.0,-29.8,-38.3,-35.3,
     $-.126,-.026,.063,.113,.208,.230,.319,.437,
     $.680,.100,.447,-.330,
     $-.360,-.359,-.443,-.420,-.383,-.375,-.329,-.294,
     $-.230,-.210,-.120,-.100,
     $-7.21,-5.40,-5.23,-5.26,-6.11,-6.64,-5.69,
     $-6.00,-6.20,-6.40,-6.60,-6.00,
     $-.380,-.363,-.378,-.386,-.370,-.453,-.550,
     $-.582,-.595,-.637,-1.02,-.840,
     $.061,.052,.052,-.012,-.013,-.024,.050,.150,
     $.130,.158,.240,.150 /
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
     $.097,.113,.145,.167,.174,.166,
     $-.048,-.038,-.040,-.021,.016,.083,
     $.127,.137,.162,.177,.179,.167,
     $-.022,-.020,-.021,-.004,.032,.094,
     $.128,.130,.154,.161,.155,.138,
     $-.040,-.038,-.039,-.025,.006,.062,
     $.087,.085,.100,.110,.104,.091,
     $-.083,-.073,-.076,-.072,-.046,.012,
     $.024,.025,.043,.053,.047,.040 /
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
     &         -1.366,-1.646,-1.917,-2.120,-2.248,-2.229 /
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
     $.245,.238,.252,.231,.198,.192,
     $.081,.077,.107,.110,.110,.141,
     $.127,.119,.133,.108,.081,.093,
     $-.046,-.020,-.009,-.005,-.006,.010,
     $.006,-.001,.014,.000,-.013,.032,
     $-.174,-.145,-.121,-.127,-.129,-.102,
     $-.097,-.113,-.087,-.084,-.069,-.006,
     $-.259,-.202,-.184,-.193,-.199,-.150,
     $-.160,-.167,-.104,-.076,-.041,-.005/
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
     $-.001,-.004,-.008,-.012,-.016,-.022,
     $-.022,-.021,-.015,-.008,-.013,-.015,
     $-.003,-.009,-.017,-.024,-.030,-.041,
     $-.045,-.040,-.016,-.002,-.010,-.019,
     $-.001,-.010,-.020,-.030,-.039,-.054,
     $-.057,-.054,-.023,-.006,-.014,-.027,
     $.000,-.010,-.022,-.034,-.047,-.060,
     $-.069,-.067,-.033,-.036,-.035,-.035,
     $.007,-.010,-.023,-.034,-.049,-.063,
     $-.081,-.079,-.060,-.058,-.062,-.059,
     $.009,-.011,-.023,-.037,-.050,-.068,
     $-.089,-.088,-.091,-.076,-.077,-.076 /
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
     $.018,.019,.018,.019,.019,.018,
     $.013,.007,.004,-.014,-.017,-.033,
     $.038,.042,.042,.042,.043,.039,
     $.030,.017,.004,-.035,-.047,-.057,
     $.056,.057,.059,.058,.058,.053,
     $.032,.012,.002,-.046,-.071,-.073,
     $.064,.077,.076,.074,.073,.057,
     $.029,.007,.012,-.034,-.065,-.041,
     $.074,.086,.093,.089,.080,.062,
     $.049,.022,.028,-.012,-.002,-.013,
     $.079,.090,.106,.106,.096,.080,
     $.068,.030,.064,.015,.011,-.001 /
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
     $-.082,-.059,-.042,-.038,-.027,-.017,
     $-.041,-.053,-.053,-.053,-.050,-.051,
     $-.066,-.043,-.038,-.027,-.023,-.016,
     $-.042,-.053,-.052,-.051,-.049,-.049,
     $-.043,-.035,-.026,-.016,-.018,-.014,
     $-.040,-.052,-.051,-.052,-.048,-.048,
     $-.042,-.037,-.031,-.026,-.017,-.012,
     $-.043,-.049,-.048,-.049,-.043,-.042,
     $-.042,-.036,-.025,-.021,-.016,-.011,
     $-.044,-.048,-.048,-.047,-.042,-.041,
     $-.020,-.028,-.013,-.014,-.011,-.010,
     $-.043,-.049,-.047,-.045,-.042,-.037,
     $-.003,-.013,-.010,-.003,-.007,-.008 /
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
     $.019,.005,.000,-.005,-.011,.008,
     $.007,.016,.014,.014,.013,.009,
     $.012,.005,.000,.004,.009,.007,
     $.013,.013,.011,.012,.011,.009,
     $.008,.005,.000,.005,.003,.005,
     $.018,.015,.015,.014,.014,.014,
     $.014,.015,.013,.011,.006,.001,
     $.015,.014,.013,.013,.012,.011,
     $.011,.010,.008,.008,.007,.003,
     $.021,.011,.010,.011,.010,.009,
     $.008,.010,.006,.005,.000,.001,
     $.023,.010,.011,.011,.011,.010,
     $.008,.010,.006,.014,.020,.000 /
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
     $.001,.017,.011,.017,.008,.016,
     $.002,-.014,-.016,-.016,-.014,-.019,
     $-.021,.002,.012,.016,.015,.011,
     $-.006,-.008,-.006,-.006,-.005,-.008,
     $-.005,.007,.004,.007,.006,.006,
     $-.011,-.011,-.010,-.009,-.008,-.006,
     $.000,.004,.007,.010,.004,.010,
     $-.015,-.015,-.014,-.012,-.011,-.008,
     $-.002,.002,.006,.012,.011,.011,
     $-.024,-.010,-.004,-.002,-.001,.003,
     $.014,.006,-.001,.004,.004,.006,
     $-.022,.002,-.003,-.005,-.003,-.001,
     $-.009,-.009,-.001,.003,-.002,.001 /
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
     $-.059,-.051,-.030,-.037,-.026,-.013,
     $-.028,-.051,-.043,-.046,-.045,-.049,
     $-.057,-.052,-.030,-.033,-.030,-.008,
     $-.037,-.041,-.038,-.040,-.040,-.038,
     $-.037,-.030,-.027,-.024,-.019,-.013,
     $-.048,-.045,-.045,-.045,-.044,-.045,
     $-.047,-.048,-.049,-.045,-.033,-.016,
     $-.043,-.044,-.041,-.041,-.040,-.038,
     $-.034,-.035,-.035,-.029,-.022,-.009,
     $-.052,-.034,-.036,-.036,-.035,-.028,
     $-.024,-.023,-.020,-.016,-.010,-.014,
     $-.062,-.034,-.027,-.028,-.027,-.027,
     $-.023,-.023,-.019,-.009,-.025,-.010 /
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


      SUBROUTINE ADC(VT,ALT,AMACH,QBAR) ! air data computer
      !This subroutine computes properties of the standard atmosphere.
      !in  : VT(m/s) ALT(m)
      !out : AMACH, QBAR(pa)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real(8) m2ft
      m2ft=3.28084
      aALT=ALT*m2ft  !ft
      aVT=VT*m2ft    !ft/s
      R0 = 0.002377d0 ! sea-level density slug/ft3
      TFAC = 1.d0 - 0.703D-5 * aALT
      T = 519.d0 * TFAC ! temperature
      IF (aALT .GE. 35000.0) T= 390.d0
      RHO = R0 * (TFAC**4.14d0)*515.3788d0 ! density (kg/m3)
      AMACH= aVT/SQRT(1.4d0*1716.3d0*T) ! Mach number
      QBAR = 0.5d0*RHO*VT*VT ! dynamic pressure (pa)
C      PS = 1715.0 * RHO * T ! static pressure
      RETURN
      END



      FUNCTION ENGINE_THRUST(THTL,ALT,RMACH)
C
C  This function computes the thrust 
C  for the F-16 model.
C
C  POW = ENGINE POWER LEVEL, PERCENT.  ( 0 <= POW <= 100. )
C  THTL= input throttle setting ( 0 <= THTL <= 1.0  
C  ALT = ALTITUDE, FT.        ( 0 <= ALT <= 50000. )
C  RMACH = MACH NUMBER.       ( 0 <= RMACH <= 1.0 )
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real(8) m2ft,lbf2N
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
      m2ft=3.28084d0
      lbf2N=4.44822d0
      H=0.0001*ALT*m2ft
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
      IF (THTL<=0.77) THEN
C
C  COMPUTE IDLE THRUST.
C
C  ALTITUDE INTERPOLATION.
        S= A(I,M)*CDH + A(I+1,M)*DH
        T= A(I,M+1)*CDH + A(I+1,M+1)*DH
C  MACH NUMBER INTERPOLATION.
        TIDL= S + (T-S)*DM
        ENGINE_THRUST= TIDL + (TMIL-TIDL)*THTL/0.77
      ELSE
C
C  COMPUTE MAX THRUST.
C
C  ALTITUDE INTERPOLATION.
        S= C(I,M)*CDH + C(I+1,M)*DH
        T= C(I,M+1)*CDH + C(I+1,M+1)*DH
C  MACH NUMBER INTERPOLATION.
        TMAX= S + (T-S)*DM
        ENGINE_THRUST= TMIL + (TMAX-TMIL)*(THTL-0.77)/(1.0-0.77)
      END IF
C
      ENGINE_THRUST=ENGINE_THRUST*lbf2N
      RETURN
      END
