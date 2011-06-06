      SUBROUTINE UACOVF ( DATA, N, RMEAN, VAR, ACVF, MLAG )
c
c    --------------------------- Official Version USHCNv2 ----------------------
c    ---------------------------------- 21 Dec 2006 ----------------------------
c
C***********************************************************************
C
C     SUBROUTINE UACOVF
C
C     FUNCTION   - UACOVF COMPUTES ESTIMATES OF THE MEAN, VARIANCE, AND
C                  AUTOCOVARIANCE FUNCTION OF A UNIVARIATE  TIME SERIES.
C
C     FORTRAN CALL - CALL UACOVF ( DATA, N, RMEAN, VAR, ACVF, MLAG )
C
C         DATA   - (REAL) AN ARRAY OF DIMENSION N. ON INPUT, THE TIME
C                  SERIES IS STORED IN DATA. ON RETURN, THE SAMPLE MEAN
C                  WILL HAVE BEEN SUBTRACTED FROM EACH ELEMENT OF DATA.
C
C         N      - (INTEGER) THE NUMBER OF OBSERVATIONS IN THE TIME
C                  SERIES.
C
C         RMEAN  - (REAL) ON OUTPUT, THE SAMPLE MEAN.
C
C         VAR    - (REAL) ON OUTPUT, THE SAMPLE VARIANCE.
C
C         ACVF   -(REAL) AN ARRAY OF DIMENSION MLAG.  ON OUTPUT THE
C                  I-TH ELEMENT  OF ACVF CONTAINS THE BIASED ESTIMATE OF
C                  THE AUTOCOVARIANCE AT LAG I.  THE ACVF AT LAG 0
C                  IS VAR.
C
C        MLAG    - (INTEGER) THE NUMBER OF LAGS OF THE AUTOCOVARIANCE
C                  THAT ARE TO BE COMPUTED.
C
C***********************************************************************
      real*8 DATA(n)
      real ACVF(n)
      DATA ZERO, IONE/0.0, 1 /
      RN = N
C
C.....COMPUTE MEAN.
C
      RMEAN = ZERO
      DO 1 I=IONE, N
         RMEAN = RMEAN + DATA(I)
 1    CONTINUE
      RMEAN = RMEAN/RN
C
C.....SUBTRACT MEAN AND COMPUTE BIASED ESTIMATE OF THE VARIANCE
C
      VAR = ZERO
      DO 2 I=IONE, N
         DATA(I) = DATA(I) - RMEAN
         VAR = VAR + DATA(I)*DATA(I)
 2    CONTINUE
      VAR = VAR/RN
C
C.....COMPUTE BIASED ESTIMATE OF AUTOCOVARIANCE FUNCTION
C
      DO 3 I = IONE, MLAG
         NMI = N - I
         SUM = ZERO
         IF (NMI) 5,5,4
 4       DO 6 J = IONE, NMI
              JPI = J+I
              SUM = SUM + DATA(J)*DATA(JPI)
 6       CONTINUE
         SUM = SUM/RN
 5       ACVF(I) = SUM
 3    CONTINUE
      RETURN
      END
