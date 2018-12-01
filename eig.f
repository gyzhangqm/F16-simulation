*     DGEEV Example Program Text
*     NAG Copyright 2005.
*     .. Parameters ..
      subroutine computeEig(AA,N)
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NB, NMAX
      PARAMETER        (NB=64,NMAX=4)
      INTEGER          LDA, LDVR, LWORK
      PARAMETER        (LDA=NMAX,LDVR=NMAX,LWORK=(2+NB)*NMAX)
*     .. Local Scalars ..
      COMPLEX *16     EIG
      INTEGER         I, INFO, J, LWKOPT, N
*     .. Local Arrays ..
      real(8)          A(LDA,NMAX),AA(LDA,NMAX)
      real(8)          DUMMY(1,1), VR(LDVR,NMAX), WI(NMAX),
     +WORK(LWORK), WR(NMAX)
*     .. External Subroutines ..
      EXTERNAL        DGEEV
*     .. Intrinsic Functions ..
      INTRINSIC       DCMPLX
*     .. Executable Statements ..
*        Compute the eigenvalues and right eigenvectors of A
*
      A=AA
      CALL DGEEV('No left vectors','Vectors (right)',N,A,LDA,WR,WI,
     +DUMMY,1,VR,LDVR,WORK,LWORK,INFO)
      LWKOPT = WORK(1)
*
      IF (INFO.EQ.0) THEN
*
*           Print solution
*
              DO 20 J = 1, N
                WRITE (NOUT,*)
                IF (WI(J).EQ.0.0D0) THEN
                  WRITE (NOUT,99999) 'Eigenvalue(', J, ') = ', WR(J)
                ELSE
                  EIG = DCMPLX(WR(J),WI(J))
                  WRITE (NOUT,99998) 'Eigenvalue(', J, ') = ', EIG
                END IF
                WRITE (NOUT,*)
                WRITE (NOUT,99997) 'Eigenvector(', J, ')'
                IF (WI(J).EQ.0.0D0) THEN
                  WRITE (NOUT,99996) (VR(I,J),I=1,N)
                ELSE IF (WI(J).GT.0.0D0) THEN
                  WRITE (NOUT,99995) (VR(I,J),VR(I,J+1),I=1,N)
                ELSE
                  WRITE (NOUT,99995) (VR(I,J-1),-VR(I,J),I=1,N)
                END IF
   20         CONTINUE
            ELSE
              WRITE (NOUT,*)
              WRITE (NOUT,99994) 'Failure in DGEEV.  INFO = ', INFO
            END IF
*
*        Print workspace information
*
      IF (LWORK.LT.LWKOPT) THEN
        WRITE (NOUT,*)
        WRITE (NOUT,99993) 'Optimum workspace required = ', LWKOPT,
     +  'Workspace provided         = ', LWORK
      END IF
*
99999 FORMAT (1X,A,I2,A,f16.9)
99998 FORMAT (1X,A,I2,A,'(',f16.9,',',f16.9,')')
99997 FORMAT (1X,A,I2,A)
99996 FORMAT (1X,f16.9)
99995 FORMAT (1X,'(',f16.9,',',f16.9,')')
99994 FORMAT (1X,A,I4)
99993 FORMAT (1X,A,I5,/1X,A,I5)
      END