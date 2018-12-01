      SUBROUTINE calAB(X,XD,U,Alongitudinal,Blongitudinal,
     &                Alateral,Blateral )
      implicit none
      integer NN,MM
      PARAMETER (NN=20, MM=10)
      real(8)::X(NN),U(MM),XD(NN)
      real(8)::V(NN),ABC(NN*NN),AAA(8,8),BBB(8,4)
      real(8)::Alongitudinal(4,4),Alateral(4,4)
      real(8)::Blongitudinal(4,2),Blateral(4,2)
      integer::NR,NC
      integer::IO(NN),JO(NN)      
      real(8)::FDX,FDU,YDX,YDU
      EXTERNAL F1,FDX,FDU,YDX,YDU
!     计算横向与纵向的A矩阵      
!            1  2     3    4   5   6   7 8 9 10        11       12 13 14 15    
!       X =  VT,ALPHA,BETA,PHI,THE,PSI,p,q,r,XN(north),XE(east),H. u, v, w
!            13 15 8  5    14 7 9 4 
!linear X =  u, w, q, THE, v, p,r,phi
      !IO是线化后选择的变量在X中的位置
      ABC(:)=0.d0; AAA(:,:)=0.d0; V(:)=0.d0
      NR=8
      NC=8
      IO(1:8)=(/13,15,8,5,14,7,9,4/)
      JO(1:8)=(/13,15,8,5,14,7,9,4/)
      V(1:15)=X(1:15)
      call JACOB (FDX,F1,X,XD,V,IO,JO,ABC,NR,NC)
      !write(*,'(8E11.4)') TRANSPOSE(reshape(ABC(1:64),(/8,8/)))
      AAA=reshape(ABC(1:64),(/8,8/))
      Alongitudinal(1:4,1:4)=AAA(1:4,1:4)
      Alateral     (1:4,1:4)=AAA(5:8,5:8)

      
      
!     计算横向与纵向的B矩阵      
!        1  2     3    4   5   6   7 8 9 10        11       12 13 14 15    
!    X = VT,ALPHA,BETA,PHI,THE,PSI,p,q,r,XN(north),XE(east),H. u, v, w
!         1   2   3   4
!    U = THTL,EL,AIL,RDR.      

      ABC(:)=0.d0; AAA(:,:)=0.d0; V(:)=0.d0
      NR=8
      NC=4
      IO(1:8)=(/13,15,8,5,14,7,9,4/)
      JO(1:4)=(/1,2,3,4/)
      V(1:4)=U(1:4)
      call JACOB (FDU,F1,X,XD,V,IO,JO,ABC,NR,NC)
      !write(*,*) 
      !write(*,'(8E11.4)') ABC(1:32)
      !write(*,*) 
      !write(*,'(4E11.4)') TRANSPOSE(reshape(ABC(1:32),(/8,4/)))
      BBB(1:8,1:4)=reshape(ABC(1:32),(/8,4/))
      Blongitudinal(1:4,1:2)=BBB(1:4,1:2)
      Blateral     (1:4,1:2)=BBB(5:8,3:4)
      
      RETURN
      END
  
      SUBROUTINE calAB1(X,XD,U,Alongitudinal,Blongitudinal,
     &                Alateral,Blateral )
      implicit none
      integer NN,MM
      PARAMETER (NN=20, MM=10)
      real(8)::X(NN),U(MM),XD(NN)
      real(8)::V(NN),ABC(NN*NN),AAA(8,8),BBB(8,4)
      real(8)::Alongitudinal(4,4),Alateral(4,4)
      real(8)::Blongitudinal(4,2),Blateral(4,2)
      integer::NR,NC
      integer::IO(NN),JO(NN)      
      real(8)::FDX,FDU,YDX,YDU
      EXTERNAL F,FDX,FDU,YDX,YDU
!     计算横向与纵向的A矩阵      
!            1  2     3    4   5   6   7 8 9 10        11       12 13 14 15    
!       X =  VT,ALPHA,BETA,PHI,THE,PSI,p,q,r,XN(north),XE(east),H. u, v, w
!            13 15 8  5    14 7 9 4 
!linear X =  u, w, q, THE, v, p,r,phi
      !IO是线化后选择的变量在X中的位置
      NR=8
      NC=8
      IO(1:8)=(/1,2,8,5,  3,7,9,4/)
      JO(1:8)=(/1,2,8,5,  3,7,9,4/)
      V(1:15)=X(1:15)
      call JACOB (FDX,F,X,XD,V,IO,JO,ABC,NR,NC)
      !write(*,'(8E11.4)') TRANSPOSE(reshape(ABC(1:64),(/8,8/)))
      AAA=reshape(ABC(1:64),(/8,8/))
      Alongitudinal(1:4,1:4)=AAA(1:4,1:4)
      Alateral     (1:4,1:4)=AAA(5:8,5:8)

      
      
!     计算横向与纵向的B矩阵      
!        1  2     3    4   5   6   7 8 9 10        11       12 13 14 15    
!    X = VT,ALPHA,BETA,PHI,THE,PSI,p,q,r,XN(north),XE(east),H. u, v, w
!         1   2   3   4
!    U = THTL,EL,AIL,RDR.      

      ABC(:)=0.d0; AAA(:,:)=0.d0; V(:)=0.d0
      NR=8
      NC=4
      IO(1:8)=(/1,2,8,5,  3,7,9,4/)
      JO(1:4)=(/1,2,3,4/)
      V(1:4)=U(1:4)
      V(2:4)=V(2:4)
      call JACOB (FDU,F,X,XD,V,IO,JO,ABC,NR,NC)
      !write(*,*) 
      !write(*,'(8E11.4)') ABC(1:32)
      !write(*,*) 
      !write(*,'(4E11.4)') TRANSPOSE(reshape(ABC(1:32),(/8,4/)))
      BBB(1:8,1:4)=reshape(ABC(1:32),(/8,4/))
      Blongitudinal(1:4,1:2)=BBB(1:4,1:2)
      Blateral     (1:4,1:2)=BBB(5:8,3:4)
      
      RETURN
      END
      SUBROUTINE calAAABBB(X,XD,U,AAA,BBB)
      implicit none
      integer NN,MM
      PARAMETER (NN=20, MM=10)
      real(8)::X(NN),U(MM),XD(NN)
      real(8)::V(NN),ABC(NN*NN),AAA(8,8),BBB(8,4)
      real(8)::Alongitudinal(4,4),Alateral(4,4)
      real(8)::Blongitudinal(4,2),Blateral(4,2)
      integer::NR,NC
      integer::IO(NN),JO(NN)      
      real(8)::FDX,FDU,YDX,YDU
      EXTERNAL F1,FDX,FDU,YDX,YDU
!     计算横向与纵向的A矩阵      
!            1  2     3    4   5   6   7 8 9 10        11       12 13 14 15    
!       X =  VT,ALPHA,BETA,PHI,THE,PSI,p,q,r,XN(north),XE(east),H. u, v, w
!            13 15 8  5    14 7 9 4 
!linear X =  u, w, q, THE, v, p,r,phi
      !IO是线化后选择的变量在X中的位置
      NR=8
      NC=8
      IO(1:8)=(/13,15,8,5,14,7,9,4/)
      JO(1:8)=(/13,15,8,5,14,7,9,4/)
      V(1:15)=X(1:15)
      call JACOB (FDX,F1,X,XD,V,IO,JO,ABC,NR,NC)
      !write(*,'(8E11.4)') TRANSPOSE(reshape(ABC(1:64),(/8,8/)))
      AAA=reshape(ABC(1:64),(/8,8/))
      Alongitudinal(1:4,1:4)=AAA(1:4,1:4)
      Alateral     (1:4,1:4)=AAA(5:8,5:8)

      
      
!     计算横向与纵向的B矩阵      
!        1  2     3    4   5   6   7 8 9 10        11       12 13 14 15    
!    X = VT,ALPHA,BETA,PHI,THE,PSI,p,q,r,XN(north),XE(east),H. u, v, w
!         1   2   3   4
!    U = THTL,EL,AIL,RDR.      

      ABC(:)=0.d0; AAA(:,:)=0.d0; V(:)=0.d0
      NR=8
      NC=4
      IO(1:8)=(/13,15,8,5,14,7,9,4/)
      JO(1:4)=(/1,2,3,4/)
      V(1:4)=U(1:4)
      V(2:4)=V(2:4)
      call JACOB (FDU,F1,X,XD,V,IO,JO,ABC,NR,NC)
      !write(*,*) 
      !write(*,'(8E11.4)') ABC(1:32)
      !write(*,*) 
      !write(*,'(4E11.4)') TRANSPOSE(reshape(ABC(1:32),(/8,4/)))
      BBB(1:8,1:4)=reshape(ABC(1:32),(/8,4/))
      Blongitudinal(1:4,1:2)=BBB(1:4,1:2)
      Blateral     (1:4,1:2)=BBB(5:8,3:4)
      
      RETURN
      END
 