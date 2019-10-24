!  =====================================================================
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!     EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
          .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMV ',INFO)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
          ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO I = 1,LENY
                      Y(I) = ZERO
                  END DO
              ELSE
                  DO I = 1,LENY
                      Y(I) = BETA*Y(I)
                  END DO
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
                  END DO
              ELSE
                  DO I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
                  END DO
              END IF
          END IF
      END IF

      IF (ALPHA.EQ.ZERO) RETURN

      IF (LSAME(TRANS,'N')) THEN
!
!        Form  y := alpha*A*x + y.
!
          JX = KX
          IF (INCY.EQ.1) THEN
              DO J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
                  END DO
                  JX = JX + INCX
              END DO
          ELSE
              DO J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
                  END DO
                  JX = JX + INCX
              END DO
          END IF

      ELSE
!
!        Form  y := alpha*A**T*x + y.
!
          JY = KY
          IF (INCX.EQ.1) THEN
              DO J = 1,N
                  TEMP = ZERO
                  DO I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
                  END DO
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
              END DO
          ELSE
              DO J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
                  END DO
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
              END DO
          END IF
      END IF
!
      RETURN
!
!     End of DGEMV .
!
      END
