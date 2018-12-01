      module var
        type :: ptr
          real(8),pointer ::p
        end type
        real(8),dimension(4) ::linearX,DeltaX,ModelX
        real(8),dimension(2) ::linearU,DeltaU
      end
      
      program main
      use var
      implicit none
      integer NN,MM
      PARAMETER (NN=20, MM=10)
      integer::NV,NX,Nstep,i,j,selection
      real(8)::COST,FDX,FDU,YDX,YDU
      real(8)::X(NN),U(MM),XD(NN),m2ft,lbf2N,RTOD,XCG
      real(8)::XD1(NN),XD2(NN)
      target::X,U
      real(8)::DT,TT,totalt
      real(8)::RADGAM,SINGAM,RR,PR,TR,PHI,CPHI,SPHI
      LOGICAL::COORD=.false., STAB=.false.
      real(8)::AN,ALAT,AX,QBAR,AMACH,Q,ALPHA,T
      real(8)::status(20,2)
      logical::Ltrim=.true.
      real(8)::V(NN),ABC(NN*NN),AAA(8,8),BBB(8,4)
      real(8)::Alongitudinal(4,4),Alateral(4,4)
      real(8)::Blongitudinal(4,2),Blateral(4,2)
      integer::NR,NC
      integer::IO(NN),JO(NN)

      COMMON/STATE/ X
      COMMON/DSTATE/ XD
      COMMON/CONTROLS/ U
      COMMON/logicalvar/Ltrim
      COMMON/PARAM/XCG
      COMMON/OUTPUT/AN,ALAT,AX,QBAR,AMACH,Q,ALPHA,T
      COMMON/CONSTRNT/RADGAM,SINGAM,RR,PR,TR,PHI,CPHI,SPHI,COORD,STAB
      EXTERNAL COST,F,FDX,FDU,YDX,YDU,F1
!                       1   2     3    4   5   6  7 8 9 10       11       12     
!    X = STATE VECTOR - VT,ALPHA,BETA,PHI,THE,PSI,P,Q,R,XN(north),XE(east),H.
!                        
C    U = INPUT VECTOR - THTL,EL,AIL,RDR.      

      m2ft=3.28084d0
      lbf2N=4.44822d0
      RTOD=180.d0/acos(-1.d0)
      XCG=0.35d0
      NV=3
      U(1:4)=(/0.2,-0.09,0.01,-0.01/) !initial guess


      !定直平飞约束
      RR=0.d0; PR=0.d0; TR=0.d0; PHI=0.d0

      call writetitle
      write(*,*) '------------------------trim-------------------------'
      !do while(.true.)
      write(*,*) 'Input Velocity(m/s),Height(m):'
      !read(*,*)  X(1),X(12)
      X(1)=130.d0; X(12)=1000.d0
      !X(1)=502.d0; X(12)=0.d0
      call TRIMMER (NV, COST)
      call printXU
      write(*,*)'推力T=',T,'N'
      !enddo


      write(*,*) '----------Numerical Linearization----------------'
      X(13)=X(1)*COS(X(2))*COS(X(3))
      X(14)=X(1)          *SIN(X(3))
      X(15)=X(1)*SIN(X(2))*COS(X(3))

      call calAB(X,XD,U,Alongitudinal,Blongitudinal,
     &                Alateral,Blateral )
      write(*,*) '纵向'
      write(*,*) 'X=[u, w, q, θ]'
      write(*,*) 'longitudinal A='
      write(*,'(5x,4f16.9)') transpose(Alongitudinal)
      write(*,*) '横向'
      write(*,*) 'X=[v, p, r, φ]'
      write(*,*) 'lateral A='
      write(*,'(5x,4f16.9)') transpose(Alateral)

      write(*,*) 'X=[u, w, q, θ] U=[THTL,EL]'
      write(*,*) '纵向longitudinal B='
      write(*,'(5x,2f16.9)') TRANSPOSE(Blongitudinal)
      write(*,*) 'X=[v, p, r, φ] U=[AIL,RDR]'
      write(*,*) '横向lateral B='
      write(*,'(5x,2f16.9)') TRANSPOSE(Blateral)
      !call printXU
      
      
      write(*,*) '----------Dynamic Responce----------------------'
      selection=1 !线性小扰动方程
      selection=2 !全量方程
333   continue
      write(*,*) '1:线性小扰动方程时域响应'
      write(*,*) '2:全量方程时域响应'
      write(*,*) '3:纵向模态激励响应'
      write(*,*) '4:横向模态激励响应'
      write(*,*) '请选择:'
      selection=1
      read(*,*)  selection
      if(selection==1) then
        goto 1
      elseif(selection==2) then
        goto 2
      elseif(selection==3) then
        goto 3
      elseif(selection==4) then
        goto 4
      else
      write(*,*) '错误数字'
      goto 333
      endif
      
1     continue
      write(*,*) '----------基于线性小扰动方程---------------------'
      selection=2
      write(*,*) 
      write(*,*) '2:升降舵单位阶跃条件状态变量的时域响应'
      write(*,*) '3:副翼单位阶跃条件状态变量的时域响应'
      write(*,*) '请选择:'
      read(*,*)  selection
      !totalt=10.d0
      write(*,*) 'Input Total Simulation time:'
      read(*,*)  totalt
      write(*,*) '仿真时间设置为',totalt,'s'
      write(*,*) '仿真时间过长可能会出现状态变量超出插值范围'
      TT=0.d0         
      DT=0.001d0 
      Nstep=int(totalt/DT)   

      select case (selection)
      case(2)

        write(*,*)'开始计算升降舵单位阶跃条件状态变量的时域响应'
!        1  2     3    4   5   6   7 8 9 10        11       12 13 14 15    
!    X = VT,ALPHA,BETA,PHI,THE,PSI,p,q,r,XN(north),XE(east),H. u, v, w
!         1   2   3   4
!    U = THTL,EL,AIL,RDR.      

        DeltaX(1)=0.d0
        DeltaX(2)=0.d0
        DeltaX(3)=0.d0
        DeltaX(4)=0.d0
        DeltaU(1)=0.d0
        DeltaU(2)=1.d0
        open(20,file='./elevator_1input_linearEq.dat')
        write(20,'(A60)')
     &  'variables=time(s),u(m/s),w(m/s),q(rad/s),THE(rad)'
        write(20,*)'ZONE T="linear"'
        do i=1,Nstep
            CALL linearF(XD1,Alongitudinal,Blongitudinal)
            do j=1,4
              DeltaX(j)=DeltaX(j)+XD1(j)*DT
            enddo
            TT=TT+DT
            linearX(1)=X(13)+DeltaX(1)
            linearX(2)=X(15)+DeltaX(2)
            linearX(3)=X(8)+DeltaX(3)
            linearX(4)=X(5)+DeltaX(4)
            linearU(1)=U(1)+DeltaU(1)
            linearU(2)=U(2)+DeltaU(2)
            write(20,'(1x,5g15.7)')TT,linearX(1:4)
          enddo
          close(20)
        case(3)
          write(*,*)'开始计算副翼单位阶跃条件状态变量的时域响应'
!        1  2     3    4   5   6   7 8 9 10        11       12 13 14 15    
!    X = VT,ALPHA,BETA,PHI,THE,PSI,p,q,r,XN(north),XE(east),H. u, v, w
!         1   2   3   4
!    U = THTL,EL,AIL,RDR.      

          DeltaX(1)=0.d0
          DeltaX(2)=0.d0
          DeltaX(3)=0.d0
          DeltaX(4)=0.d0
          DeltaU(1)=1.d0
          DeltaU(2)=0.d0
          open(20,file='./ailerons_1input_linearEq.dat')
          write(20,'(A60)')
     &    'variables=time(s),v(m/s),p(rad/s),r(rad/s),PHI(rad)'
          write(20,*)'ZONE T="linear"'
          do i=1,Nstep
              CALL linearF(XD1,Alateral,Blateral)
              do j=1,4
                DeltaX(j)=DeltaX(j)+XD1(j)*DT
              enddo
              TT=TT+DT
              linearX(1)=X(14)+DeltaX(1)
              linearX(2)=X(7)+DeltaX(2)
              linearX(3)=X(9)+DeltaX(3)
              linearX(4)=X(4)+DeltaX(4)
              write(20,'(1x,5g15.7)')TT,linearX(1:4)
            enddo
            close(20)
          end select
        
      write(*,*) '仿真结束，请根据生成文件作图'
      write(*,*) '按任意键结束'
      read(*,*)
      stop      
      
      
2     continue
      write(*,*) '----------基于全量运动方程-----------------------'
      selection=2
      write(*,*) 
      write(*,*) '1:零输入下状态变量的时域响应'
      write(*,*) '2:升降舵单位阶跃条件状态变量的时域响应'
      write(*,*) '3:副翼单位阶跃条件状态变量的时域响应'
      write(*,*) '请选择:'
      read(*,*)  selection
      !totalt=10.d0
      write(*,*) 'Input Total Simulation time:'
      read(*,*)  totalt
      write(*,*) '仿真时间设置为',totalt,'s'
      write(*,*) '仿真时间过长可能会出现状态变量超出插值范围'
      DT=0.001d0 
      Nstep=int(totalt/DT)   

      select case (selection)
      case (1)
        write(*,*)'开始计算零输入下状态变量的时域响应'
        TT=0.d0           
        NX=12           
        open(10,file='./zero_input.dat')
        write(10,'(A120)')'variables=time(s),u(m/s),v(m/s),w(m/s),
     &PHI(rad),THE(rad),PSI(rad),p(rad/s),q(rad/s),r(rad/s),
     &xE(m),yE(m),H(m)'
        write(10,*)'ZONE T="non-linear"'
        do i=1,Nstep
          call RK4(F,DT,X,NX)
          TT=TT+DT
          write(10,'(1x,13g15.7)')TT,X(1)*cos(X(2))*cos(X(3)),
     &       X(1)*sin(X(3)),X(1)*sin(X(2))*cos(X(3)),
     &       X(4:9),X(10:12)
        enddo
        close(10)

      case (2)
        write(*,*)'开始计算升降舵单位阶跃条件状态变量的时域响应'
        U(2)=U(2)+1.d0
        TT=0.d0         
        NX=12           
        open(20,file='./elevator_1input.dat')
        write(20,'(A120)')'variables=time(s),u(m/s),v(m/s),w(m/s),
     &  PHI(rad),THE(rad),PSI(rad),p(rad/s),q(rad/s),r(rad/s),
     &  xE(m),yE(m),H(m)'
        write(20,*)'ZONE T="non-linear"'
        do i=1,Nstep
          !X(3)=0.0         
          call RK4(F,DT,X,NX)
          TT=TT+DT
          write(20,'(1x,13g15.7)')TT,X(1)*cos(X(2))*cos(X(3)),
     &       X(1)*sin(X(3)),X(1)*sin(X(2))*cos(X(3)),
     &       X(4:9),X(10:12)
        enddo
        close(20)
      case (3)
        write(*,*)'开始计算副翼单位阶跃条件下状态变量的时域响应'
        U(3)=U(3)+1.d0
        TT=0.d0     
        NX=12   
        open(30,file='./ailerons_1input.dat')
        write(30,'(A120)')'variables=time(s),u(m/s),v(m/s),w(m/s),
     &  PHI(rad),THE(rad),PSI(rad),p(rad/s),q(rad/s),r(rad/s),
     &  xE(m),yE(m),H(m)'
        write(30,*)'ZONE T="non-linear"'
        do i=1,Nstep
          !X(3)=0.0         
          call RK4(F,DT,X,NX)
          TT=TT+DT
          write(30,'(1x,13g15.7)')TT,X(1)*cos(X(2))*cos(X(3)),
     &       X(1)*sin(X(3)),X(1)*sin(X(2))*cos(X(3)),
     &       X(4:9),X(10:12)
        enddo
        close(30)
      end select

      write(*,*) '仿真结束，请根据生成文件作图'
      write(*,*) '按任意键结束'
      read(*,*)
      stop

3     continue
      write(*,*) '----------纵向模态激励响应---------------------'
      
      call computeEig(Alongitudinal,4)!计算特征值及特征向量
      write(*,*)'请输入一个特征值对应的特征向量的实部'
      read(*,*)ModelX(1:4)
      !totalt=10.d0
      write(*,*) 'Input Total Simulation time:'
      read(*,*)  totalt
      write(*,*) '仿真时间设置为',totalt,'s'
      write(*,*)'开始计算升降舵单位阶跃条件状态变量的时域响应'
      TT=0.d0         
      DT=0.001d0 
      Nstep=int(totalt/DT)   
!        1  2     3    4   5   6   7 8 9 10        11       12 13 14 15    
!    X = VT,ALPHA,BETA,PHI,THE,PSI,p,q,r,XN(north),XE(east),H. u, v, w
!         1   2   3   4
!    U = THTL,EL,AIL,RDR.      

        DeltaX(1:4)=ModelX(1:4)
        DeltaU(:)=0.d0
        open(20,file='./zongxiangmotaijili.dat')
        write(20,'(A120)')
     &  'variables=time(s),u(m/s),w(m/s),q(rad/s),THE(rad)'
        do i=1,Nstep
            call linearRK4(DT,Alongitudinal,Blongitudinal)
            !CALL linearF(XD1,Alongitudinal,Blongitudinal)
            !do j=1,4
            !  DeltaX(j)=DeltaX(j)+XD1(j)*DT
            !enddo
            linearX(1)=X(13)+DeltaX(1)
            linearX(2)=X(15)+DeltaX(2)
            linearX(3)=X(8)+DeltaX(3)
            linearX(4)=X(5)+DeltaX(4)
            TT=TT+DT
            write(20,'(1x,5g15.7)')TT,linearX(1:4)
          enddo
          close(20)     
      
      write(*,*) '仿真结束，请根据生成文件作图'
      write(*,*) '按任意键结束'
      read(*,*)
      stop
     
      
4     continue
      write(*,*) '----------横向模态激励响应---------------------'
      call computeEig(Alateral,4)!计算特征值及特征向量
      write(*,*)'请输入一个特征值对应的特征向量的实部'
      read(*,*)ModelX(1:4)
      !totalt=10.d0
      write(*,*) 'Input Total Simulation time:'
      read(*,*)  totalt
      write(*,*) '仿真时间设置为',totalt,'s'
      TT=0.d0         
      DT=0.001d0 
      Nstep=int(totalt/DT)   
!        1  2     3    4   5   6   7 8 9 10        11       12 13 14 15    
!    X = VT,ALPHA,BETA,PHI,THE,PSI,p,q,r,XN(north),XE(east),H. u, v, w
!         1   2   3   4
!    U = THTL,EL,AIL,RDR.      

        DeltaX(1:4)=ModelX(1:4)
        DeltaU(:)=0.d0
        open(20,file='./hengxiangmotaijili.dat')
        write(20,'(A120)')
     &    'variables=time(s),v(m/s),p(rad/s),r(rad/s),PHI(rad)'
        do i=1,Nstep
            call linearRK4(DT,Alateral,Blateral)
            linearX(1)=X(14)+DeltaX(1)
            linearX(2)=X(7)+DeltaX(2)
            linearX(3)=X(9)+DeltaX(3)
            linearX(4)=X(4)+DeltaX(4)
            TT=TT+DT
            write(20,'(1x,5g15.7)')TT,linearX(1:4)
          enddo
          close(20)     
      
      write(*,*) '仿真结束，请根据生成文件作图'
      write(*,*) '按任意键结束'
      read(*,*)
      stop
      end





      subroutine main_test_SMPLX
      implicit none
      integer NN,MM
      PARAMETER (NN=20, MM=10)

      integer NV,NC
      real(8) COST,bananafunction
      real(8) X(NN),U(MM),m2ft,S(2),DS(2),SIGMA,F0,FFIN
      EXTERNAL COST,bananafunction
!                       1   2     3    4   5   6  7 8 9 10       11       12     
!    X = STATE VECTOR - VT,ALPHA,BETA,PHI,THE,PSI,P,Q,R,XN(orth),XE(east),H.
C    U = INPUT VECTOR - THTL,EL,AIL,RDR.      
      COMMON/STATE/ X
      COMMON/ CONTROLS/ U

      m2ft=3.28084
      NV=3
      X=0.0
      U=0.0
      X(1) =400.0        !Vt
      X(12)=0.0       


      !call TRIMMER (NV, COST)
      !test 
      nv=2        
      S=(/2.0,0.0/)   
      DS=(/0.5,0.5/)
      SIGMA=-1    
      NC=10000
      CALL SMPLX(bananafunction,NV,S,DS,SIGMA,NC,F0,FFIN)
      WRITE(*,*)F0,FFin,S
      end


      FUNCTION bananafunction (S) !
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION S(*)
      bananafunction = (S(1)-1)**2+5*(S(1)**2-S(2))**2
      RETURN
      END

      subroutine printXU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer NN,MM
      PARAMETER (NN=20, MM=10)
      real(8)::X(NN),U(MM)
      COMMON/STATE/ X
      COMMON/CONTROLS/ U      
      RTOD=180.d0/acos(-1.d0)
      WRITE(*,'(/1X,A)')'  Alpha°       Throttle    Elevator°
     &    Ailerons°    Rudder° '
      WRITE(*,'(1X,5(1PE10.2,3X))')
     &      X(2)*RTOD,U(1),U(2),X(3)*RTOD, U(4)*RTOD
        
      write(*,*) "X:"
      write(*,*) "  u =",X(1)*cos(X(2))*cos(X(3))
      write(*,*) "  v =",X(1)*sin(X(3))
      write(*,*) "  w =",X(1)*sin(X(2))*cos(X(3))
      write(*,*) " φ =",X(4)*RTOD,'°'
      write(*,*) " θ =",X(5)*RTOD,'°'
      write(*,*) " ψ =",X(6)*RTOD,'°'
      write(*,*) "  p =",X(7),'rad/s'
      write(*,*) "  q =",X(8),'rad/s'
      write(*,*) "  r =",X(9),'rad/s'
      write(*,*) "  xE=",X(10),'m'
      write(*,*) "  yE=",X(11),'m'
      write(*,*) "  zE=",X(12),'m'
      RETURN
      END

	subroutine writetitle
	character*10 bar10
	character*20 bar20
	bar10='----------'
	bar20='--------------------'
	write(6,'(4a20)')bar20,bar20,bar20,bar20
	write(6,*)'        Non-linear F16 Simulation            '
	write(6,*)'             开发者：张启明'
	write(6,'(4a20)')bar20,bar20,bar20,bar20
	write(6,*)'  重要提示：'
	write(6,*)'::本程序根据Aircraft Control and Simulation课本改写'
	write(6,*)'::正确性已与
     & http://www.aem.umn.edu/~balas/darpa_sec/SEC.Software.html
     & 上的程序进行了对比验证 '
	!write(6,*)'::老师给定的气动数据是作用在X=0.35c处的'
	!write(6,*)'::F16的质心是在X=0.30c处'
	write(6,*)'::当计算气动系数时，副翼偏角和方向舵偏角也需要无量纲化'
	write(6,*)'::分别使用其舵偏角最大值进行无量纲化'
	write(6,*)'::本程序暂时只对全量运动方程进行时域响应计算'
	write(6,*)'::本程序对线化方程使用数值微分的方法计算'
	write(6,*)'::其它事项可与开发者联系.'
	write(6,'(4a20)')bar20,bar20,bar20,bar20
	write(6,*)

	end      