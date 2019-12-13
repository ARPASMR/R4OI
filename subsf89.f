C SUBROUTINES
C xzinv
C tdxzinv
C a5inv
C FUNCTIONS
C det5
C det4
C*******************************************************************************
      SUBROUTINE XZINV(KK, YO, SLON, SHGT, AVLON, DZBF
     &               , ZINV, YINV, ALPA, GAMA, ALPB, GAMB, RESMIN)
C-------------------------------------------------------------------------------
C height and longitude
C YBk= GAMA*(SHGT(K)-ZINV)+YINV  above
C YBk= GAMB*(SHGT(K)-ZINV)+YINV  below
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER KK
      REAL*4 SLON(KK), SHGT(KK), YO(KK)
C-----------------------------------------------------------------------
      REAL*8 AVLON
      REAL*8 TINV, DXA, DZA, DXB, DZB
C
      REAL*8 AVHGT
      REAL*8 ABK, ABT
      REAL*8 AX, AZ, AX2, AZ2, AZX, AXT, AZT
      REAL*8 BX, BZ, BX2, BZ2, BZX, BXT, BZT
C-----------------------------------------------------------------------
      REAL*8 AA(5,5), ANV(5,5), DET
      REAL*8 YB, RES, RESMIN
C-----------------------------------------------------------------------
      REAL*8 ZINV, YINV, ALPA, GAMA, ALPB, GAMB
     &, RES1, ALP1, GAM1, AVT
      REAL*8 DZBF
      REAL*8 BFABOV, BFBELO, ZABOV, ZBELO
C-----------------------------------------------------------------------
C      INTEGER I, J
      INTEGER K, INV, KHA, KHB
     &, ICUTMX, ICUTMN, ICTD, ICUT
      INTEGER KHADEF, KHBDEF
C      CHARACTER ANS*1
C-----------------------------------------------------------------------
      KHADEF= 0
      KHBDEF= 0
C-----------------------------------------------------------------------
      ABK= KK
      ABT= 0.
      AVLON= 0.
      AVHGT= 0.
      DO K=1,KK
       ABT= ABT+ YO(K)
       AVLON= AVLON+ SLON(K)
       AVHGT= AVHGT+ SHGT(K)
      END DO
      AVLON= AVLON /ABK
      AVHGT= AVHGT /ABK
      AVT= ABT/ ABK
C      WRITE (6,*) 'xzinv: ', KK, ABK, AVLON, AVHGT
C       WRITE (6,*) ' ok?'
C       READ (5,'(A1)') ANS
C      WRITE (6,*) ' '
C-------------------------------------------------------------------------------
C prima faccio il conto senza inversione :
C-------------------------------------------------------------------------------
       WRITE (6,*) ' costruisco yb senza inversione'
       AX2= 0.
       AZ2= 0.
       AZX= 0.
       AXT= 0.
       AZT= 0.
C-----------------------------------------------------------------------
       DO K=1,KK
        AX2= AX2+ (SLON(K)-AVLON)**2
        AZ2= AZ2+ (SHGT(K)-AVHGT)**2
        AZX= AZX+  (SHGT(K)-AVHGT) *(SLON(K)-AVLON)
        AXT= AXT+ (YO(K)-AVT) *(SLON(K) -AVLON)
        AZT= AZT+ (YO(K)-AVT) *(SHGT(K) -AVHGT)
       END DO
       DET= AX2*AZ2-AZX*AZX
       WRITE (6,*) 'xzinv : s.i. ax2 az2 azx axt azt det:'
     &, AX2, AZ2, AZX, AXT, AZT, DET
C-----------------------------------------------------------------------
       ALP1= (AZ2*AXT -AZX*AZT) /DET
       GAM1= (AX2*AZT -AZX*AXT) /DET
C-----------------------------------------------------------------------
       RES1= 0.
       DO K=1,KK
        YB= AVT +ALP1*(SLON(K)-AVLON) +GAM1*(SHGT(K)-AVHGT)
        RES1= RES1 +(YB-YO(K))*(YB-YO(K))
       END DO
       RES1= SQRT(RES1/KK)
       WRITE (6,*) 'xzinv : senza inversione res1:',RES1 
       WRITE (6,*) 'xzinv : s.i. avT avZ Tx Tz '
     &, AVT, AVHGT, ALP1, GAM1
C       WRITE (6,*) ' ok?'
C       READ (5,'(A1)') ANS
C-------------------------------------------------------------------------------
      ICUTMX= 2000
      ICUTMN= 0
      ICTD= -10
      RESMIN= 100000.
C-------------------------------------------------------------------------------
      DO ICUT=ICUTMX,ICUTMN,ICTD
C       WRITE (6,*) ' icut', ICUT
C-----------------------------------------------------------------------
C A
C-----------------------------------------------------------------------
       AX= 0.
       AZ= 0.
       AX2= 0.
       AZ2= 0.
       AZX= 0.
       AXT= 0.
       AZT= 0.
       KHA= 0
C-----------------------------------------------------------------------
C B
C-----------------------------------------------------------------------
       BX= 0.
       BZ= 0.
       BX2= 0.
       BZ2= 0.
       BZX= 0.
       BXT= 0.
       BZT= 0.
       KHB= 0
C-----------------------------------------------------------------------
       DO K=1,KK
        IF (SHGT(K).GT.ICUT) THEN
         KHA= KHA +1
         AX =  AX+ (SLON(K)-AVLON)
         AZ =  AZ+  (SHGT(K)-ICUT)
         AX2= AX2+ (SLON(K)-AVLON)**2
         AZ2= AZ2+  (SHGT(K)-ICUT)**2
         AZX= AZX+   (SHGT(K)-ICUT) *(SLON(K)-AVLON)
         AXT= AXT+ YO(K) *(SLON(K) -AVLON)
         AZT= AZT+ YO(K) *(SHGT(K) -ICUT)
        ELSE
         KHB= KHB +1
         BX =  BX+ (SLON(K)-AVLON)
         BZ =  BZ+  (SHGT(K)-ICUT)
         BX2= BX2+ (SLON(K)-AVLON)**2
         BZ2= BZ2+  (SHGT(K)-ICUT)**2
         BZX= BZX+   (SHGT(K)-ICUT) *(SLON(K)-AVLON)
         BXT= BXT+ YO(K) *(SLON(K) -AVLON)
         BZT= BZT+ YO(K) *(SHGT(K) -ICUT)
        END IF
       END DO
C       WRITE (6,*) 'kk kha khb', KK,KHA,KHB
C       WRITE (6,*) ' ok?'
C       READ (5,'(A1)') ANS
       IF (KHA.LE.6) CYCLE
       IF (KHB.LE.6) EXIT
C-----------------------------------------------------------------------
       AA(1,1)= ABK 
       AA(1,2)= AX
       AA(1,3)= AZ
       AA(1,4)= BX
       AA(1,5)= BZ
C
       AA(2,1)= AA(1,2)
       AA(2,2)= AX2
       AA(2,3)= AZX
       AA(2,4)= 0.
       AA(2,5)= 0.
C
       AA(3,1)= AA(1,3)
       AA(3,2)= AA(2,3)
       AA(3,3)= AZ2
       AA(3,4)= 0.
       AA(3,5)= 0.
C
       AA(4,1)= AA(1,4)
       AA(4,2)= AA(2,4)
       AA(4,3)= AA(3,4)
       AA(4,4)= BX2
       AA(4,5)= BZX
C
       AA(5,1)= AA(1,5)
       AA(5,2)= AA(2,5)
       AA(5,3)= AA(3,5)
       AA(5,4)= AA(4,5)
       AA(5,5)= BZ2
C-----------------------------------------------------------------------
C       WRITE (6,*) 'matrice'
C       DO I=1,5
C        WRITE (6,'(5G11.3)') (AA(I,J),J=1,5)
C       END DO
       INV= 0
       CALL A5INV(AA, DET, ANV)
       IF (DET.EQ.0) THEN
        WRITE (6,*) 'xyzinv : matrix is not invertible'
        WRITE (6,*) 'xzinv : kha khb icut ', KHA, KHB, ICUT
        STOP
       END IF
C       WRITE (6,*) 'inversa'
C       DO I=1,5
C        WRITE (6,'(5G11.3)') (ANV(I,J),J=1,5)
C       END DO
C-----------------------------------------------------------------------
       TINV= ANV(1,1)*ABT +ANV(1,2)*AXT +ANV(1,3)*AZT
     &                    +ANV(1,4)*BXT +ANV(1,5)*BZT
        DXA= ANV(2,1)*ABT +ANV(2,2)*AXT +ANV(2,3)*AZT
     &                    +ANV(2,4)*BXT +ANV(2,5)*BZT
        DZA= ANV(3,1)*ABT +ANV(3,2)*AXT +ANV(3,3)*AZT
     &                    +ANV(3,4)*BXT +ANV(3,5)*BZT
        DXB= ANV(4,1)*ABT +ANV(4,2)*AXT +ANV(4,3)*AZT
     &                    +ANV(4,4)*BXT +ANV(4,5)*BZT
        DZB= ANV(5,1)*ABT +ANV(5,2)*AXT +ANV(5,3)*AZT
     &                    +ANV(5,4)*BXT +ANV(5,5)*BZT
C-----------------------------------------------------------------------
C calcolo residui sui dati
C-----------------------------------------------------------------------
       RES= 0.
C raccordo:
       ZABOV= ICUT +DZBF
       ZBELO= ICUT -DZBF
C
       DO K=1,KK
        IF (SHGT(K).GT.ZABOV) THEN
         YB= TINV +DXA *(SLON(K) -AVLON) +DZA *(SHGT(K) -ICUT)
        ELSE IF (SHGT(K).LT.ZBELO) THEN
         YB= TINV +DXB *(SLON(K) -AVLON) +DZB *(SHGT(K) -ICUT)
        ELSE
         BFABOV= TINV +DXA *(SLON(K) -AVLON) +DZA *(ZABOV -ICUT)
         BFBELO= TINV +DXB *(SLON(K) -AVLON) +DZB *(ZBELO -ICUT)
         YB= (BFABOV *(SHGT(K) -ZBELO) +BFBELO *(ZABOV -SHGT(K)))
     &      /(ZABOV -ZBELO)
        END IF
        RES= RES +(YB-YO(K))*(YB-YO(K))
       END DO
       RES= SQRT(RES/KK)
C-----------------------------------------------------------------------
C       IF (RESMIN.GE.RES) THEN
C-----------------------------------------------------------------------
       IF (ABS(DXA).LE.2. .AND. ABS(DXB).LE.2. .AND.
     &     RESMIN.GE.RES) THEN
C-----------------------------------------------------------------------
C anche condizione su gradiente verticale "above"?
C -.011 deg/m= -11. deg/km :adiabatica secca
C          .AND. DZA.GE.-0.011 .AND. DZA.LE.0.0
C-----------------------------------------------------------------------
        RESMIN= RES
        ZINV= ICUT
        YINV= TINV
        ALPA= DXA
        GAMA= DZA
        ALPB= DXB
        GAMB= DZB
        KHADEF= KHA
        KHBDEF= KHB
       END IF
C       WRITE (6,*) 'icut kha khb yinv res resmin zinv '
C     &            , ICUT,KHA,KHB,TINV,RES,RESMIN,ZINV
C       WRITE (6,'(A15,3F15.8)') 'above Tx Tz ', DXA,DZA
C       WRITE (6,'(A15,3F15.8)') 'below Tx Tz ', DXB,DZB
C       WRITE (6,*) ' '
C       WRITE (6,*) ' ok?'
C       READ (5,'(A1)') ANS
C-----------------------------------------------------------------------
      END DO
C-------------------------------------------------------------------------------
C
      WRITE(6,*) 'numero di stazioni sopra l''inversione (KHADEF)'
     &, KHADEF
      WRITE(6,*) 'numero di stazioni sotto l''inversione (KHBDEF)'
     &, KHBDEF
C
C      IF (RES1.LE.RESMIN .OR. GAMA.GT.0) THEN
      IF (RES1.LE.RESMIN) THEN
       WRITE (6,*) 'xzinv: MEGLIO SENZA'
       ZINV= AVHGT
       YINV= AVT
       ALPA= ALP1
       ALPB= ALP1
       GAMA= GAM1
       GAMB= GAM1
      END IF
C      WRITE (6,*) ' end xzinv'
C      WRITE (6,*) ' ok?'
C      READ (5,'(A1)') ANS
C-------------------------------------------------------------------------------
      RETURN
      END
C*******************************************************************************
      SUBROUTINE TDXZINV(KK, YO, STA, SLON, SHGT, AVLON, DZBF
     &                 , ZINV, YINV, ALPA, GAMA, ALPB, GAMB, RESMIN)
C-------------------------------------------------------------------------------
C height and longitude
C YBk= GAMA*(SHGT(K)-ZINV)+YINV  above
C YBk= GAMB*(SHGT(K)-ZINV)+YINV  below
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER KK
      REAL*4 SLON(KK), SHGT(KK), YO(KK), STA(KK)
C-----------------------------------------------------------------------
      REAL*8 AVLON
      REAL*8 TINV, DXA, DZA, DXB, DZB
C
      REAL*8 AVHGT
      REAL*8 ABK, ABT
      REAL*8 AX, AZ, AX2, AZ2, AZX, AXT, AZT
      REAL*8 BX, BZ, BX2, BZ2, BZX, BXT, BZT
C-----------------------------------------------------------------------
      REAL*8 AA(5,5), ANV(5,5), DET
      REAL*8 YB, RES, RESMIN
C-----------------------------------------------------------------------
      REAL*8 ZINV, YINV, ALPA, GAMA, ALPB, GAMB
     &, RES1, ALP1, GAM1, AVT
      REAL*8 DZBF
      REAL*8 BFABOV, BFBELO, ZABOV, ZBELO
C maximum oversaturated Td value:
      REAL*4 TDMAX
C-----------------------------------------------------------------------
C      INTEGER I, J
      INTEGER K, INV, KHA, KHB
     &, ICUTMX, ICUTMN, ICTD, ICUT
      INTEGER KHADEF, KHBDEF
C      CHARACTER ANS*1
C-----------------------------------------------------------------------
      KHADEF= 0
      KHBDEF= 0
C-----------------------------------------------------------------------
      ABK= KK
      ABT= 0.
      AVLON= 0.
      AVHGT= 0.
      DO K=1,KK
       ABT= ABT+ YO(K)
       AVLON= AVLON+ SLON(K)
       AVHGT= AVHGT+ SHGT(K)
      END DO
      AVLON= AVLON /ABK
      AVHGT= AVHGT /ABK
      AVT= ABT/ ABK
C      WRITE (6,*) 'xzinv: ', KK, ABK, AVLON, AVHGT
C       WRITE (6,*) ' ok?'
C       READ (5,'(A1)') ANS
C      WRITE (6,*) ' '
C-------------------------------------------------------------------------------
C prima faccio il conto senza inversione :
C-------------------------------------------------------------------------------
       WRITE (6,*) ' costruisco yb senza inversione'
       AX2= 0.
       AZ2= 0.
       AZX= 0.
       AXT= 0.
       AZT= 0.
C-----------------------------------------------------------------------
       DO K=1,KK
        AX2= AX2+ (SLON(K)-AVLON)**2
        AZ2= AZ2+ (SHGT(K)-AVHGT)**2
        AZX= AZX+  (SHGT(K)-AVHGT) *(SLON(K)-AVLON)
        AXT= AXT+ (YO(K)-AVT) *(SLON(K) -AVLON)
        AZT= AZT+ (YO(K)-AVT) *(SHGT(K) -AVHGT)
       END DO
       DET= AX2*AZ2-AZX*AZX
       WRITE (6,*) 'xzinv : s.i. ax2 az2 azx axt azt det:'
     &, AX2, AZ2, AZX, AXT, AZT, DET
C-----------------------------------------------------------------------
       ALP1= (AZ2*AXT -AZX*AZT) /DET
       GAM1= (AX2*AZT -AZX*AXT) /DET
C-----------------------------------------------------------------------
       RES1= 0.
       DO K=1,KK
        YB= AVT +ALP1*(SLON(K)-AVLON) +GAM1*(SHGT(K)-AVHGT)
        RES1= RES1 +(YB-YO(K))*(YB-YO(K))
       END DO
       RES1= SQRT(RES1/KK)
       WRITE (6,*) 'xzinv : senza inversione res1:',RES1 
       WRITE (6,*) 'xzinv : s.i. avT avZ Tx Tz '
     &, AVT, AVHGT, ALP1, GAM1
C       WRITE (6,*) ' ok?'
C       READ (5,'(A1)') ANS
C-------------------------------------------------------------------------------
      ICUTMX= 2000
      ICUTMN= 0
      ICTD= -10
      RESMIN= 100000.
C-------------------------------------------------------------------------------
      DO ICUT=ICUTMX,ICUTMN,ICTD
C       WRITE (6,*) ' icut', ICUT
C-----------------------------------------------------------------------
C A
C-----------------------------------------------------------------------
       AX= 0.
       AZ= 0.
       AX2= 0.
       AZ2= 0.
       AZX= 0.
       AXT= 0.
       AZT= 0.
       KHA= 0
C-----------------------------------------------------------------------
C B
C-----------------------------------------------------------------------
       BX= 0.
       BZ= 0.
       BX2= 0.
       BZ2= 0.
       BZX= 0.
       BXT= 0.
       BZT= 0.
       KHB= 0
C-----------------------------------------------------------------------
       DO K=1,KK
        IF (SHGT(K).GT.ICUT) THEN
         KHA= KHA +1
         AX =  AX+ (SLON(K)-AVLON)
         AZ =  AZ+  (SHGT(K)-ICUT)
         AX2= AX2+ (SLON(K)-AVLON)**2
         AZ2= AZ2+  (SHGT(K)-ICUT)**2
         AZX= AZX+   (SHGT(K)-ICUT) *(SLON(K)-AVLON)
         AXT= AXT+ YO(K) *(SLON(K) -AVLON)
         AZT= AZT+ YO(K) *(SHGT(K) -ICUT)
        ELSE
         KHB= KHB +1
         BX =  BX+ (SLON(K)-AVLON)
         BZ =  BZ+  (SHGT(K)-ICUT)
         BX2= BX2+ (SLON(K)-AVLON)**2
         BZ2= BZ2+  (SHGT(K)-ICUT)**2
         BZX= BZX+   (SHGT(K)-ICUT) *(SLON(K)-AVLON)
         BXT= BXT+ YO(K) *(SLON(K) -AVLON)
         BZT= BZT+ YO(K) *(SHGT(K) -ICUT)
        END IF
       END DO
C       WRITE (6,*) 'kk kha khb', KK,KHA,KHB
C       WRITE (6,*) ' ok?'
C       READ (5,'(A1)') ANS
       IF (KHA.LE.6) CYCLE
       IF (KHB.LE.6) EXIT
C-----------------------------------------------------------------------
       AA(1,1)= ABK 
       AA(1,2)= AX
       AA(1,3)= AZ
       AA(1,4)= BX
       AA(1,5)= BZ
C
       AA(2,1)= AA(1,2)
       AA(2,2)= AX2
       AA(2,3)= AZX
       AA(2,4)= 0.
       AA(2,5)= 0.
C
       AA(3,1)= AA(1,3)
       AA(3,2)= AA(2,3)
       AA(3,3)= AZ2
       AA(3,4)= 0.
       AA(3,5)= 0.
C
       AA(4,1)= AA(1,4)
       AA(4,2)= AA(2,4)
       AA(4,3)= AA(3,4)
       AA(4,4)= BX2
       AA(4,5)= BZX
C
       AA(5,1)= AA(1,5)
       AA(5,2)= AA(2,5)
       AA(5,3)= AA(3,5)
       AA(5,4)= AA(4,5)
       AA(5,5)= BZ2
C-----------------------------------------------------------------------
C       WRITE (6,*) 'matrice'
C       DO I=1,5
C        WRITE (6,'(5G11.3)') (AA(I,J),J=1,5)
C       END DO
       INV= 0
       CALL A5INV(AA, DET, ANV)
       IF (DET.EQ.0) THEN
        WRITE (6,*) 'xyzinv : matrix is not invertible'
        WRITE (6,*) 'xzinv : kha khb icut ', KHA, KHB, ICUT
        STOP
       END IF
C       WRITE (6,*) 'inversa'
C       DO I=1,5
C        WRITE (6,'(5G11.3)') (ANV(I,J),J=1,5)
C       END DO
C-----------------------------------------------------------------------
       TINV= ANV(1,1)*ABT +ANV(1,2)*AXT +ANV(1,3)*AZT
     &                    +ANV(1,4)*BXT +ANV(1,5)*BZT
        DXA= ANV(2,1)*ABT +ANV(2,2)*AXT +ANV(2,3)*AZT
     &                    +ANV(2,4)*BXT +ANV(2,5)*BZT
        DZA= ANV(3,1)*ABT +ANV(3,2)*AXT +ANV(3,3)*AZT
     &                    +ANV(3,4)*BXT +ANV(3,5)*BZT
        DXB= ANV(4,1)*ABT +ANV(4,2)*AXT +ANV(4,3)*AZT
     &                    +ANV(4,4)*BXT +ANV(4,5)*BZT
        DZB= ANV(5,1)*ABT +ANV(5,2)*AXT +ANV(5,3)*AZT
     &                    +ANV(5,4)*BXT +ANV(5,5)*BZT
C-----------------------------------------------------------------------
C calcolo residui sui dati
C-----------------------------------------------------------------------
       RES= 0.
C raccordo:
       ZABOV= ICUT +DZBF
       ZBELO= ICUT -DZBF
C
       DO K=1,KK
        IF (SHGT(K).GT.ZABOV) THEN
         YB= TINV +DXA *(SLON(K) -AVLON) +DZA *(SHGT(K) -ICUT)
        ELSE IF (SHGT(K).LT.ZBELO) THEN
         YB= TINV +DXB *(SLON(K) -AVLON) +DZB *(SHGT(K) -ICUT)
        ELSE
         BFABOV= TINV +DXA *(SLON(K) -AVLON) +DZA *(ZABOV -ICUT)
         BFBELO= TINV +DXB *(SLON(K) -AVLON) +DZB *(ZBELO -ICUT)
         YB= (BFABOV *(SHGT(K) -ZBELO) +BFBELO *(ZABOV -SHGT(K)))
     &      /(ZABOV -ZBELO)
        END IF
C YBSAT
C check for saturation against T analysis values
C no oversaturated values in the background field:
        TDMAX= STA(K)
        IF (YB.GT.TDMAX) YB= TDMAX
C
        RES= RES +(YB-YO(K))*(YB-YO(K))
       END DO
       RES= SQRT(RES/KK)
C-----------------------------------------------------------------------
C       IF (RESMIN.GE.RES) THEN
C-----------------------------------------------------------------------
       IF (ABS(DXA).LE.2. .AND. ABS(DXB).LE.2. .AND.
     &     RESMIN.GE.RES) THEN
C-----------------------------------------------------------------------
C anche condizione su gradiente verticale "above"?
C -.011 deg/m= -11. deg/km :adiabatica secca
C          .AND. DZA.GE.-0.011 .AND. DZA.LE.0.0
C-----------------------------------------------------------------------
        RESMIN= RES
        ZINV= ICUT
        YINV= TINV
        ALPA= DXA
        GAMA= DZA
        ALPB= DXB
        GAMB= DZB
        KHADEF= KHA
        KHBDEF= KHB
       END IF
C       WRITE (6,*) 'icut kha khb yinv res resmin zinv '
C     &            , ICUT,KHA,KHB,TINV,RES,RESMIN,ZINV
C       WRITE (6,'(A15,3F15.8)') 'above Tx Tz ', DXA,DZA
C       WRITE (6,'(A15,3F15.8)') 'below Tx Tz ', DXB,DZB
C       WRITE (6,*) ' '
C       WRITE (6,*) ' ok?'
C       READ (5,'(A1)') ANS
C-----------------------------------------------------------------------
      END DO
C-------------------------------------------------------------------------------
C
      WRITE(6,*) 'numero di stazioni sopra l''inversione (KHADEF)'
     &, KHADEF
      WRITE(6,*) 'numero di stazioni sotto l''inversione (KHBDEF)'
     &, KHBDEF
C
C      IF (RES1.LE.RESMIN .OR. GAMA.GT.0) THEN
      IF (RES1.LE.RESMIN) THEN
       WRITE (6,*) 'xzinv: MEGLIO SENZA'
       ZINV= AVHGT
       YINV= AVT
       ALPA= ALP1
       ALPB= ALP1
       GAMA= GAM1
       GAMB= GAM1
      END IF
C      WRITE (6,*) ' end xzinv'
C      WRITE (6,*) ' ok?'
C      READ (5,'(A1)') ANS
C-------------------------------------------------------------------------------
      RETURN
      END
C*******************************************************************************
      SUBROUTINE A5INV(A,DET,ANV)
C-------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=5)
      REAL*8 A(N,N), ANV(N,N), B(N-1,N-1)
      REAL*8 DET, DET4, DET5
      INTEGER I, J, II, JJ, IB, JB, ISIGN
C-----------------------------------------------------------------------
      DET= 0.
      DET= DET5(A)
C      WRITE (6,*) 'A5inv det=',DET
C      IF (DET.EQ.0) RETURN 'non invertibile'
      IF (DET.EQ.0) RETURN
C-----------------------------------------------------------------------
      DO I=1,N
       DO J=1,N
        ISIGN= (-1)**(I+J)
C-----------------------------------------------------------------------
C costruzione matrice ausiliaria B
C-----------------------------------------------------------------------
        IB= 0
        DO II=1,N
         IF (II.EQ.I) CYCLE
         IB= IB +1
         JB= 0
         DO JJ=1,N
          IF (JJ.EQ.J) CYCLE
          JB= JB +1
          B(IB,JB)= A(II,JJ)
         END DO
        END DO
C        WRITE (6,*) 'i j sign',I,J,ISIGN,' ok?'
        ANV(J,I)= ISIGN* DET4(B) /DET
C-----------------------------------------------------------------------
       END DO
      END DO
C-------------------------------------------------------------------------------
      RETURN
      END
C*******************************************************************************
      REAL*8 FUNCTION DET5(A)
C-------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=5)
      REAL*8 A(N,N), B(N-1,N-1), DET4
      INTEGER ISIGN, K, I, J, IB, JB
C-----------------------------------------------------------------------
      DET5= 0.
      ISIGN= -1
      DO K=1,N
       JB= 0
       DO J=1,N
        IF (J.NE.K) THEN
         JB= JB +1
         DO I=2,N
          IB= I -1
          B(IB,JB) = A(I,J)
         END DO
        END IF
       END DO
       ISIGN= -ISIGN
       DET5= DET5 +ISIGN*A(1,K)*DET4(B)
      END DO
C-------------------------------------------------------------------------------
      RETURN
      END
C*******************************************************************************
      REAL*8 FUNCTION DET4(A)
C-------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 A(4,4)
      REAL*8 DET11, DET12, DET13, DET14
C-------------------------------------------------------------------------------
      DET11=  A(1,1) *( A(2,2)* (A(3,3)*A(4,4) -A(3,4)*A(4,3))
     &                 +A(2,3)* (A(3,4)*A(4,2) -A(3,2)*A(4,4))
     &                 +A(2,4)* (A(3,2)*A(4,3) -A(3,3)*A(4,2)) )
C-----------------------------------------------------------------------
      DET12= -A(1,2) *( A(2,1)* (A(3,3)*A(4,4) -A(3,4)*A(4,3))
     &                 +A(2,3)* (A(3,4)*A(4,1) -A(3,1)*A(4,4))
     &                 +A(2,4)* (A(3,1)*A(4,3) -A(3,3)*A(4,1)) )
C-----------------------------------------------------------------------
      DET13=  A(1,3) *( A(2,1)* (A(3,2)*A(4,4) -A(3,4)*A(4,2))
     &                 +A(2,2)* (A(3,4)*A(4,1) -A(3,1)*A(4,4))
     &                 +A(2,4)* (A(3,1)*A(4,2) -A(3,2)*A(4,1)) )
C-----------------------------------------------------------------------
      DET14= -A(1,4) *( A(2,1)* (A(3,2)*A(4,3) -A(3,3)*A(4,2))
     &                 +A(2,2)* (A(3,3)*A(4,1) -A(3,1)*A(4,3))
     &                 +A(2,3)* (A(3,1)*A(4,2) -A(3,2)*A(4,1)) )
C-----------------------------------------------------------------------
      DET4= DET11 +DET12 +DET13 +DET14
C-------------------------------------------------------------------------------
      RETURN
      END
C*******************************************************************************
