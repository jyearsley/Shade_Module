PROGRAM RIPARIAN
  IMPLICIT NONE
!
  CHARACTER (LEN = 1)   :: EQUALS
  CHARACTER (LEN = 5)   :: CCLOCK
  CHARACTER (LEN = 8)   :: DDATE,DUMMY1
  INTEGER               :: I,II,MM,NDAY,XI
  REAL                  :: JDAY,LDAY,LAT,LONGIT,PHI0
  REAL                  :: EXTIN,ISW,S_FCTR,TREE_HT,W_BANK,W_BUFFER,W_STREAM
  REAL,DIMENSION(12)    :: EXTIN_COEFF
  INTEGER,DIMENSION(12) :: MONTH = [31,59,90,120,151,181,212,243,273,304,334,365]
!
  OPEN (UNIT = 10,FILE = 'Test_ISW.txt',STATUS = 'OLD')
  OPEN (UNIT = 20,FILE = 'Stream_Params',STATUS = 'OLD')
    READ(20,*) DUMMY1,EQUALS,LONGIT
    write(*,*) 'LONGIT - ',LONGIT
    READ(20,*) DUMMY1,EQUALS,LAT
    write(*,*) 'LONGIT - ',LAT
    READ(20,*) DUMMY1,EQUALS,PHI0
    write(*,*) 'LONGIT - ',PHI0
    READ(20,*) DUMMY1,EQUALS,TREE_HT
    write(*,*) 'LONGIT - ',TREE_HT
    READ(20,*) DUMMY1,EQUALS,W_BANK
    write(*,*) 'LONGIT - ',W_BANK
    READ(20,*) DUMMY1,EQUALS,W_BUFFER
    write(*,*) 'LONGIT - ',W_BUFFER
    READ(20,*) DUMMY1,EQUALS,W_STREAM
    write(*,*) 'LONGIT - ',W_STREAM
    READ(20,*) DUMMY1,EQUALS,(EXTIN_COEFF(I),I=1,12)
    write(*,*) 'LONGIT - ',EXTIN_COEFF
    
    
  OPEN (UNIT = 30,FILE = 'SHADE.Out',STATUS = 'UNKNOWN') 
!
!
  MM   = 1
  NDAY = 1
  JDAY = 1.0
!
  DO WHILE (JDAY .LE. 365.0)
    EXTIN = EXTIN_COEFF(MM)
    DO II = 1,8
      XI = II
      LDAY = JDAY + (3.0*(XI-0.5))/24.  
      READ (10,*) DDATE,CCLOCK,ISW
      CALL SHADING (LDAY,LONGIT,LAT,PHI0,TREE_HT,W_BANK,W_BUFFER,W_STREAM,EXTIN,S_FCTR)
      write(30,*) JDAY,LDAY,EXTIN,S_FCTR
    END DO
    JDAY = JDAY + 1
    IF (JDAY .GT. MONTH(MM)) MM = MM + 1
  END DO
END PROGRAM RIPARIAN
!***********************************************************************************************************************************
!**                                                S U B R O U T I N E   S H A D I N G                                            **
!***********************************************************************************************************************************

SUBROUTINE SHADING (JDAY,LONGIT,LAT,PHI0,TREE_HT,W_BANK,W_BUFFER,W_STREAM,EXTIN,SHADE_FCTR)
  IMPLICIT NONE
  CHARACTER(1)    :: BANK
  REAL            :: A0COS_FCTR,A0SIN_FCTR,A0TAN_FCTR
  REAL            :: AZ_FCTR,AZCOS_FCTR,AZSIN_FCTR
  REAL            :: EQTNEW,LOCAL,STANDARD,HOUR,TAUD,SINAL,LAT_RAD,A02,AZ00,A0,AX
  REAL            :: A00,A00_DEG,DECL,HH,LAT,LONGIT,PHI0
  REAL            :: ANG1,ANG2,TOPOANG,AZT           ! Inactive until topographic shading is included-JRY
  REAL            :: TREE_HT,JDAY,W_STREAM,W_BANK,W_BANK_AZ,W_BUFFER,W_NOSHADE,W_SHADE
  REAL            :: EXTIN,PL1,PL2,PL_AVG,SHADE_FCTR,STLEN
  REAL            :: DX1,DY1,DX2,DY2
  INTEGER         :: CASE_NO
  INTEGER         :: IDAY,J
  LOGICAL         :: SHADE = .TRUE.
  REAL,PARAMETER  :: DEG_RAD = 3.1415927/180.,RAD_DEG = 180./3.1415927,PI = 3.1415927
!
    PHI0     = DEG_RAD*PHI0
    LAT_RAD  = DEG_RAD*LAT
!
! Calculate solar altitude, declination, and local hour angle when short-wave solar radiation is provided as input

    LOCAL    =  LONGIT
    STANDARD =  15.0*INT(LONGIT/15.0)
    HOUR     = (JDAY-INT(JDAY))*24.0
    IDAY     =  JDAY-((INT(JDAY/365))*365)
    IDAY     =  IDAY+INT(INT(JDAY/365)/4)
    TAUD     = (2*PI*(IDAY-1))/365
    EQTNEW   =  0.170*SIN(4*PI*(IDAY-80)/373)-0.129*SIN(2*PI*(IDAY-8)/355)
    HH       =  0.261799*(HOUR-(LOCAL-STANDARD)*0.0666667+EQTNEW-12.0)
    DECL     =  0.006918-0.399912*COS(TAUD)+0.070257*SIN(TAUD)-0.006758*COS(2*TAUD)+0.000907*SIN(2*TAUD)-0.002697*COS(3*TAUD)      &
                +0.001480*SIN(3*TAUD)
    SINAL    =  SIN(LAT_RAD)*SIN(DECL)+COS(LAT_RAD)*COS(DECL)*COS(HH)
    A00      =  ASIN(SINAL)
    A00_DEG  =  RAD_DEG*A00
! If the sun is below the horizon, set SHADE(I) to 0

  IF (A00_DEG < 0.0) THEN
    CASE_NO = 0
  END IF

!** Calculate solar azimuth angle

    A02 = A00
    AX  = (SIN(DECL)*COS(LAT_RAD)-COS(DECL)*COS(HH)*SIN(LAT_RAD))/COS(A02)
 write(50,*) JDAY,IDAY,HH,DECL*RAD_DEG,A00_DEG,AX*RAD_DEG
    IF (AX >  1.0) AX =  1.0
    IF (AX < -1.0) AX = -1.0
    AZT = ACOS(AX)
    IF (HH < 0.0) THEN
     AZ00 = AZT
    ELSE
     AZ00 = 2.0*PI-AZT
    END IF
    A0 = A02

!** Interpolate the topographic shade angle
!
!    DO J=1,IANG-1
!      IF (AZ00 > ANG(J) .AND. AZ00 <= ANG(J+1)) THEN
!        ANG1    =  AZ00-ANG(J)
!        ANG2    = (TOPO(I,J+1)-TOPO(I,J))/GAMA                 ! SW 10/17/05
!        TOPOANG =  TOPO(I,J)+ANG2*ANG1
!      END IF
!    END DO
!    IF (AZ00 > ANG(IANG) .AND. AZ00 <= 2*PI) THEN
!      ANG1    =  AZ00-ANG(IANG)
!      ANG2    = (TOPO(I,1)-TOPO(I,IANG))/GAMA                  ! SW 10/17/05
!      TOPOANG =  TOPO(I,IANG)+ANG2*ANG1
!    END IF

!** Complete topographic shading if solar altitude less than topo angle
!
!   IF (A0 <= TOPOANG) THEN
!      SFACT = 0.90
!      GO TO 100
!    END IF

!** Bank with the controlling vegetation

    IF (PHI0 > 0.0 .AND. PHI0 <= PI) THEN
      IF (AZ00 > PHI0     .AND. AZ00 <= PHI0+PI) BANK = 'L'
      IF (AZ00 > 0.0         .AND. AZ00 <= PHI0)    BANK = 'R'
      IF (AZ00 > PHI0+PI  .AND. AZ00 <  2.0*PI)     BANK = 'R'
    ELSE IF (PHI0 > PI .AND. PHI0 <= 2.0*PI) THEN
      IF (AZ00 >= PHI0    .AND. AZ00 < 2.0*PI)      BANK = 'L'
      IF (AZ00 >= 0.0        .AND. AZ00 < PHI0-PI)  BANK = 'L'
      IF (AZ00 >= PHI0-PI .AND. AZ00 < PHI0)     BANK = 'R'
    END IF

!   Estimate shade for left or right bank. For the time being, assume
!  
    Select Case (BANK)
!
!   Select riparian shading data from the left bank
!
      Case ('L')
!
!       Initialize the shading properties for the left bank
! 
!
!   Select riparian shading data from the right bank
!
      Case ('R')
!
!       Initialize the shading properties for the right bank
! 
    End Select
!
!   Calculate some shading properties
!
!    W_Shade = (Tree_Ht*SIN(PHI0-AZ00)/tan(A00)


!** Distance from vegetation to water edge on line parallel to azimuth

!    W_BANK_AZ = W_BANK/ABS(SIN(PHI0-AZ00))
!    IF (STLEN <= W_BANK_AZ) THEN
!      SHADE_FCTR = 0.0
!     GO TO 100
!    END IF

!** No vegetative shading if azimuth angle is oriented parallel to stream

    IF (AZ00 == PHI0 .OR. AZ00 == PHI0+PI .OR. AZ00+PI == PHI0) THEN
      CASE_NO = 1
    END IF
!
!   Calculate some important parameters for determining riparian shading
!
    AZ_FCTR = ABS(PHI0-AZ00)
    AZCOS_FCTR = COS(AZ_FCTR)
    AZSIN_FCTR = SIN(AZ_FCTR)
    A0COS_FCTR = COS(A00)
    A0SIN_FCTR = SIN(A00)
    A0TAN_FCTR = TAN(A00)
!
    DX1 = (TREE_HT/TAN(A00))-(W_STREAM+W_BANK)
    DY1 = TREE_HT-(W_STREAM+W_BANK)*TAN(A00)
!
    DX2 = (TREE_HT/TAN(A00)) - W_BANK
    DY2 = TREE_HT - (W_BANK)*TAN(A00)
!
    PL1 = 0.0
!
!   Partial and no shade - Cases 1, 2 & 3
!
    IF (DX1 .LT. 0.0) THEN 
!  
      IF (DX2 .LT. 0.0) CASE_NO = 1
!
      IF (DX2 .GE. 0.0 .AND. DX2 .LT. W_BUFFER-0.1) CASE_NO = 2
!
      IF (DX2 .GE. W_BUFFER-1) CASE_NO = 3
!
   END IF    
!
!   Full shade - Cases 4, 5 & 6
!
   IF (DX1 .GE. 0.0 .AND. DX1 .LT. W_BUFFER-0.1) THEN
!
     IF (DX2 .LT. W_BUFFER-0.1) CASE_NO = 4
!
     IF (DX2 .GE. W_BUFFER-0.1) CASE_NO = 5
   END IF
!
   IF (DX1 .GE. W_BUFFER-0.1) CASE_NO = 6
!
   SELECT CASE(CASE_NO)
!
     CASE (1)
!
       SHADE_FCTR = 0.0
       write(50,*) CASE_NO,A00*RAD_DEG,SHADE_FCTR
!
     CASE (2)
!
       DY2 = TREE_HT - (W_STREAM + W_BANK)*TAN(A00)
       W_SHADE = (TREE_HT*AZSIN_FCTR/A0TAN_FCTR) - W_BANK
       W_SHADE = ABS(W_SHADE)
       W_NOSHADE = W_STREAM - W_SHADE
       PL2 = DY2/(A0SIN_FCTR*AZSIN_FCTR)
       PL2 = ABS(PL2)
       PL_AVG = 0.5*(PL1 + PL2)
!
       SHADE_FCTR = (W_SHADE + EXP(-EXTIN*PL_AVG)*W_NOSHADE)/W_STREAM
       write(50,*) CASE_NO,AZ00*RAD_DEG,PL_AVG,W_SHADE,W_NOSHADE,SHADE_FCTR
!       
     CASE (3)
!
       DY2 = (W_STREAM + W_BANK)*TAN(A00)
       W_SHADE = (TREE_HT*AZSIN_FCTR/A0TAN_FCTR) - W_BANK
       W_SHADE = ABS(W_SHADE)
       W_NOSHADE = W_STREAM - W_SHADE
       PL2 = DY2/(A0SIN_FCTR*AZSIN_FCTR)
       PL2 = ABS(PL2)
   
       PL_AVG = 0.5*(PL1 + PL2)
!
       SHADE_FCTR = (W_SHADE + EXP(-EXTIN*PL_AVG)*W_NOSHADE)/W_STREAM
       write(50,*) CASE_NO,AZ00*RAD_DEG,PL_AVG,W_SHADE,W_NOSHADE,SHADE_FCTR
!
     CASE (4)
!
       DY1 = TREE_HT - (W_STREAM + W_BANK) * A0TAN_FCTR
       PL1 = DY1/(A0SIN_FCTR*AZSIN_FCTR)
       PL1 = ABS(PL1)
!
       DY2 = TREE_HT - W_BANK * A0TAN_FCTR
       PL2 = DY2/(A0SIN_FCTR*AZSIN_FCTR)
       PL2 = ABS(PL2)
       PL_AVG = 0.5 * (PL1 + PL2)
!
       SHADE_FCTR = EXP(-EXTIN*PL_AVG)       
       write(50,*) CASE_NO,AZ00*RAD_DEG,PL_AVG,SHADE_FCTR
!
     CASE (5)
!
       DY1 = TREE_HT - (W_STREAM + W_BANK) * A0TAN_FCTR
       PL1 = DY1/(A0SIN_FCTR*AZSIN_FCTR)
       PL1 = ABS(PL1)
!
       DY2 = W_BUFFER * A0TAN_FCTR
       PL2 = DY2/(A0SIN_FCTR*AZSIN_FCTR)
       PL2 = ABS(PL2)
       PL_AVG = 0.5 * (PL1 + PL2)
!
       SHADE_FCTR = EXP(-EXTIN*PL_AVG) 
       write(50,*) CASE_NO,AZ00*RAD_DEG,PL_AVG,SHADE_FCTR

!
     CASE (6)
!
       DY1 = W_BUFFER * A0TAN_FCTR
       DY2 = DY1
       PL_AVG = DY1/(A0SIN_FCTR*AZSIN_FCTR)
       PL_AVG = ABS(PL_AVG)
!
       SHADE_FCTR = EXP(-EXTIN*PL_AVG) 
       write(50,*) CASE_NO,AZ00*RAD_DEG,PL_AVG,SHADE_FCTR
!
   END SELECT       
!     

!** Path factor consider solar altitude, declination of sun and solar altitude
!
!    P_Factor = SIN(PHI0(NSEG)-AZ00)/TAN(A0)

!    SN    = MIN (HT*ABS (SIN (ABS (PHI0(I)-AZ00)))/TAN (A0)-WDBANK,BI(KT,I))
!    SFACT = SRED*SN/BI(KT,I)
100 CONTINUE
!    SHADE_FCTR = MAX (0.0,1-SHADE_FCTR)
  RETURN
END SUBROUTINE SHADING
