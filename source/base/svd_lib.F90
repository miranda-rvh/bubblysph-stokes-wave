module svd_lib
  use kind_parameters
  implicit none

!  private

contains


!! ------------------------------------------------------------------------------------------------
  subroutine svd_solve(AMATRX,NSIZE,BVECTR)

!C     DESCRIPTION
!C     -----------
!C     COMPUTES THE SINGULAR VALUE DECOMPOSITION OF A MATRIX
!C
!C     IN THE FORM A = U.W.V^T
!C     
!C     INPUT:
!C       MATRIX AMATRX(NSIZE,NSIZE) OF PHYSICAL DIMENSIONS nsize,nsize
!C
!C     OUTPUT:
!C       MATRIX U(NSIZE,NSIZE) STORED IN AMATRX
!C       MATRIX W = DIAG(NSIZE) STORED AS WVECTR OF PHYSICAL DIMENSION nsize
!C       MATRIX V(NSIZE,NSIZE) (NOT THE TRANSPOSE) STORED IN VMATRX(nsize,nsize)
!C
!C     BASED ON NUMERICAL RECIPES SUBROUTINE SVDCMP
!C
!C     *************************************************************************

!C     PARAMETERS
!C     ==========
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      INTEGER NMAX
      PARAMETER(NMAX=500)

!C     ARGUMENTS
!C     =========
      DOUBLE PRECISION, intent(in) :: AMATRX(nsize,nsize)
      DOUBLE PRECISION VMATRX(nsize,nsize),UMATRX(nsize,nsize)
      DOUBLE PRECISION WVECTR(nsize)
      INTEGER NSIZE


!C     ARGUMENTS
!C     =========
      DOUBLE PRECISION, intent(inout) :: BVECTR(nsize)



!C     EXTERNAL FUNCTION
!C     =================
!      DOUBLE PRECISION PYTHAG
!      EXTERNAL PYTHAG


!C     LOCAL DATA
!C     ==========
      DOUBLE PRECISION TMP(NMAX)
      DOUBLE PRECISION RV1(NMAX)
      DOUBLE PRECISION ANORM,CVAR,FVAR,GVAR,HVAR,SVAR,SSCALE
      DOUBLE PRECISION XVAR,YVAR,ZVAR
      INTEGER IC,ITS,JC,JJ,KC,LC,NM


!!    Set umatrx to amatrx
     umatrx = amatrx

!C     BEGIN
!C     =====
!C     HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM
      GVAR = ZERO
      SSCALE = ZERO
      ANORM = ZERO
      DO IC = 1,NSIZE
        LC = IC + 1
        RV1(IC) = SSCALE*GVAR
        GVAR = ZERO
        SVAR = ZERO
        SSCALE = ZERO
        IF(IC.LE.NSIZE)THEN
          DO KC = IC,NSIZE
            SSCALE = SSCALE + ABS(uMATRX(KC,IC))
          ENDDO
          IF(SSCALE.NE.ZERO)THEN
            DO KC = IC,NSIZE
              uMATRX(KC,IC) = uMATRX(KC,IC)/SSCALE
              SVAR = SVAR + uMATRX(KC,IC)*uMATRX(KC,IC)
            ENDDO
            FVAR = uMATRX(IC,IC)
            GVAR = -SIGN(SQRT(SVAR),FVAR)
            HVAR = FVAR*GVAR - SVAR
            uMATRX(IC,IC) = FVAR - GVAR
            DO JC = LC,NSIZE
              SVAR = ZERO
              DO KC = IC,NSIZE
                SVAR = SVAR + uMATRX(KC,IC)*uMATRX(KC,JC)
              ENDDO
              FVAR = SVAR/HVAR
              DO KC = IC,NSIZE
                uMATRX(KC,JC) = uMATRX(KC,JC) + FVAR*uMATRX(KC,IC)
              ENDDO
            ENDDO
            DO KC = IC,NSIZE
              uMATRX(KC,IC) = SSCALE*uMATRX(KC,IC)
            ENDDO
          ENDIF
        ENDIF
        WVECTR(IC) = SSCALE*GVAR
        GVAR = ZERO
        SVAR = ZERO
        SSCALE = ZERO
        IF(IC.LE.NSIZE)THEN
          DO KC = LC,NSIZE
            SSCALE = SSCALE + ABS(uMATRX(IC,KC))
          ENDDO
          IF(SSCALE.NE.ZERO)THEN
            DO KC = LC,NSIZE
              uMATRX(IC,KC) = uMATRX(IC,KC)/SSCALE
              SVAR = SVAR + uMATRX(IC,KC)*uMATRX(IC,KC)
            ENDDO
            FVAR = uMATRX(IC,LC)
            GVAR = -SIGN(SQRT(SVAR),FVAR)
            HVAR = FVAR*GVAR - SVAR
            uMATRX(IC,LC) = FVAR - GVAR
            DO KC = LC,NSIZE
              RV1(KC)=uMATRX(IC,KC)/HVAR
            ENDDO
            DO JC = LC,NSIZE
              SVAR = ZERO
              DO KC = LC,NSIZE
                SVAR = SVAR + uMATRX(JC,KC)*uMATRX(IC,KC)
              ENDDO
              DO KC = LC,NSIZE
                uMATRX(JC,KC) = uMATRX(JC,KC) + SVAR*RV1(KC)
              ENDDO
            ENDDO
            DO KC = LC,NSIZE
              uMATRX(IC,KC)=SSCALE*uMATRX(IC,KC)
            ENDDO
          ENDIF
        ENDIF
        ANORM = MAX(ANORM,(ABS(WVECTR(IC))+ABS(RV1(IC))))
      ENDDO

!C     ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS
      DO IC = NSIZE,1,-1
        IF(IC.LT.NSIZE)THEN
          IF(GVAR.NE.ZERO)THEN
!C           DOUBLE DIVISION TO AVOID POSSIBLE UNDERFLOW
            DO JC = LC,NSIZE
              VMATRX(JC,IC) = (uMATRX(IC,JC)/uMATRX(IC,LC))/GVAR
            ENDDO
            DO JC = LC,NSIZE
              SVAR = ZERO
              DO KC = LC,NSIZE
                SVAR = SVAR + uMATRX(IC,KC)*VMATRX(KC,JC)
              ENDDO
              DO KC = LC,NSIZE
                VMATRX(KC,JC) = VMATRX(KC,JC) + SVAR*VMATRX(KC,IC)
              ENDDO
            ENDDO
          ENDIF
          DO JC = LC,NSIZE
            VMATRX(IC,JC) = ZERO
            VMATRX(JC,IC) = ZERO
          ENDDO
        ENDIF
        VMATRX(IC,IC) = ONE
        GVAR = RV1(IC)
        LC = IC
      ENDDO

!C     ACCUMULATION OF LEFT-HAND TRANSFORMATIONS
      DO IC = NSIZE,1,-1
        LC = IC + 1
        GVAR = WVECTR(IC)
        DO JC = LC,NSIZE
          uMATRX(IC,JC) = ZERO
        ENDDO
        IF(GVAR.NE.ZERO)THEN
          GVAR = ONE/GVAR
          DO JC=LC,NSIZE
            SVAR = ZERO
            DO KC = LC,NSIZE
              SVAR = SVAR + uMATRX(KC,IC)*uMATRX(KC,JC)
            ENDDO
            FVAR = (SVAR/uMATRX(IC,IC))*GVAR
            DO KC = IC,NSIZE
              uMATRX(KC,JC) = uMATRX(KC,JC) + FVAR*uMATRX(KC,IC)
            ENDDO
          ENDDO
          DO JC = IC,NSIZE
            uMATRX(JC,IC) = uMATRX(JC,IC)*GVAR
          ENDDO
        ELSE
          DO JC= IC,NSIZE
            uMATRX(JC,IC) = ZERO
          ENDDO
        ENDIF
        uMATRX(IC,IC) = uMATRX(IC,IC) + ONE
      ENDDO

!C     DIAGONALIZATION OF THE BIDIAGONAL FORM
!C     LOOP OVER SINGULAR VALUES
      NM = 1
      DO KC = NSIZE,1,-1
!C       LOOP OVER ALLOWED ITERATIONS
        DO ITS = 1,30
          DO LC = KC,1,-1
!C           TEST FOR SPLITTING
            NM = LC-1
!C           NOTE THAT RV1(1) IS ALWAYS ZERO
            IF((ABS(RV1(LC))+ANORM).EQ.ANORM) GOTO 2000
            IF((ABS(WVECTR(NM))+ANORM).EQ.ANORM) GOTO 1000
          ENDDO

!C         CANCELLATION OF RV1(LC) IF LC > 1
1000      CVAR = ZERO
          SVAR = ONE
          DO IC = LC,KC
            FVAR = SVAR*RV1(IC)
            RV1(IC) = CVAR*RV1(IC)
            IF((ABS(FVAR)+ANORM).EQ.ANORM) GOTO 2000
            GVAR = WVECTR(IC)
            HVAR = PYTHAG(FVAR,GVAR)
            WVECTR(IC) = HVAR
            HVAR = ONE/HVAR
            CVAR =  (GVAR*HVAR)
            SVAR = -(FVAR*HVAR)
            DO JC = 1,NSIZE
              YVAR = uMATRX(JC,NM)
              ZVAR = uMATRX(JC,IC)
              uMATRX(JC,NM) =  (YVAR*CVAR)+(ZVAR*SVAR)
              uMATRX(JC,IC) = -(YVAR*SVAR)+(ZVAR*CVAR)
            ENDDO
          ENDDO

2000      ZVAR = WVECTR(KC)
          IF(LC.EQ.KC)THEN
!C           CONVERGENCE
            IF(ZVAR.LT.ZERO)THEN
!C             SINGULAR VALUE IS MADE NONNEGATIVE
              WVECTR(KC) = -ZVAR
              DO JC = 1,NSIZE
                VMATRX(JC,KC) = -VMATRX(JC,KC)
              ENDDO
            ENDIF
            GOTO 3000
          ENDIF
          IF(ITS.EQ.30)THEN
            WRITE(6,*)'SVDCMP: no convergence'
            WRITE(6,'(4I7)')NSIZE
            DO IC = 1,NSIZE
              WRITE(6,'(19(1PE15.7))')(uMATRX(IC,JC),JC=1,NSIZE)
            ENDDO
            WRITE(6,*)
            DO IC = 1,NSIZE
              WRITE(6,'(19(1PE15.7))')WVECTR(IC)
            ENDDO
            WRITE(6,*)
            DO IC = 1,NSIZE
              WRITE(6,'(19(1PE15.7))')(VMATRX(IC,JC),JC=1,NSIZE)
            ENDDO
            WRITE(6,*)
            STOP
          ENDIF

!C         SHIFT FROM BOTTOM 2-BY-2 MINOR
          XVAR = WVECTR(LC)
          NM = KC-1
          YVAR = WVECTR(NM)
          GVAR = RV1(NM)
          HVAR = RV1(KC)
          FVAR = ((YVAR-ZVAR)*(YVAR+ZVAR) +  (GVAR-HVAR)*(GVAR+HVAR))/(TWO*HVAR*YVAR)
          GVAR = PYTHAG(FVAR,ONE)
          FVAR = ((XVAR-ZVAR)*(XVAR+ZVAR) +  HVAR*((YVAR/(FVAR+SIGN(GVAR,FVAR)))-HVAR))/XVAR
!C         NEXT QR TRANSFORMATION
          CVAR = ONE
          SVAR = ONE
          DO JC = LC,NM
            IC = JC+1
            GVAR = RV1(IC)
            YVAR = WVECTR(IC)
            HVAR = SVAR*GVAR
            GVAR = CVAR*GVAR
            ZVAR = PYTHAG(FVAR,HVAR)
            RV1(JC) = ZVAR
            CVAR = FVAR/ZVAR
            SVAR = HVAR/ZVAR
            FVAR =  (XVAR*CVAR)+(GVAR*SVAR)
            GVAR = -(XVAR*SVAR)+(GVAR*CVAR)
            HVAR = YVAR*SVAR
            YVAR = YVAR*CVAR
            DO JJ = 1,NSIZE
              XVAR = VMATRX(JJ,JC)
              ZVAR = VMATRX(JJ,IC)
              VMATRX(JJ,JC) =  (XVAR*CVAR)+(ZVAR*SVAR)
              VMATRX(JJ,IC) = -(XVAR*SVAR)+(ZVAR*CVAR)
            ENDDO
            ZVAR = PYTHAG(FVAR,HVAR)
            WVECTR(JC) = ZVAR
!C           ROTATION CAN BE ARBITRARY IF ZVAR = 0
            IF(ZVAR.NE.ZERO)THEN
              ZVAR = ONE/ZVAR
              CVAR = FVAR*ZVAR
              SVAR = HVAR*ZVAR
            ENDIF
            FVAR =  (CVAR*GVAR)+(SVAR*YVAR)
            XVAR = -(SVAR*GVAR)+(CVAR*YVAR)
            DO JJ = 1,NSIZE
              YVAR = uMATRX(JJ,JC)
              ZVAR = uMATRX(JJ,IC)
              uMATRX(JJ,JC) =  (YVAR*CVAR)+(ZVAR*SVAR)
              uMATRX(JJ,IC) = -(YVAR*SVAR)+(ZVAR*CVAR)
            ENDDO
          ENDDO
          RV1(LC) = ZERO
          RV1(KC) = FVAR
          WVECTR(KC) = XVAR
        ENDDO

3000    CONTINUE

      ENDDO


!     *************************************************************************
      !! Now use the SVD to solve the linear system

!     BEGIN
!     =====
!     CALCULATE U^T.B
      DO JC = 1,NSIZE
        SVAR = ZERO
!       NONZERO RESULT ONLY IF W(J) IS NONZERO
        IF(WVECTR(JC).NE.ZERO)THEN
          DO IC = 1,NSIZE
            SVAR = SVAR + UMATRX(IC,JC)*BVECTR(IC)
          ENDDO
!         THIS IS THE DIVIDE BY W(J) .
          SVAR = SVAR/WVECTR(JC)
        ENDIF
        TMP(JC) = SVAR
      ENDDO

!     MATRIX MULTIPLY BY V TO GET THE ANSWER
      DO JC = 1,NSIZE
        SVAR = ZERO
        DO JJ = 1,NSIZE
          SVAR = SVAR + VMATRX(JC,JJ)*TMP(JJ)
        ENDDO
        BVECTR(JC) = SVAR
      ENDDO

      RETURN
   end subroutine svd_solve


   FUNCTION PYTHAG(AVALUE,BVALUE)

!C     *************************************************************************
!C
!C     PYTHAG
!C     ======
!C     DESCRIPTION
!C     -----------
!C     COMPUTES (A^2 + B^2)^1/2 WITHOUT DESTRUCTIVE UNDERFLOW OR OVERFLOW
!C
!C     BASED ON NUMERICAL RECIPES FUNCTION PYTHAG
!C
!C     *************************************************************************

!C     PARAMETERS
!C     ==========
      DOUBLE PRECISION ZERO,ONE
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)


!C     FUNCTION VALUE
!C     ==============
      DOUBLE PRECISION PYTHAG


!C     ARGUMENTS
!C     =========
      DOUBLE PRECISION AVALUE,BVALUE


!C     LOCAL DATA
!C     ==========
      DOUBLE PRECISION ABSA,ABSB


!C     BEGIN
!C     =====
      ABSA = ABS(AVALUE)
      ABSB = ABS(BVALUE)
      IF(ABSA.GT.ABSB)THEN
        PYTHAG = ABSA*SQRT(ONE+(ABSB/ABSA)**2)
      ELSE
        IF(ABSB.EQ.ZERO)THEN
          PYTHAG = ZERO
        ELSE
          PYTHAG = ABSB*SQRT(ONE+(ABSA/ABSB)**2)
        ENDIF
      ENDIF


      RETURN
      END function pythag


!! ------------------------------------------------------------------------------------------------
end module svd_lib
