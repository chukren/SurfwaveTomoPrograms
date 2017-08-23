	subroutine s_saito
C PROGRAM FOR DISPERSION PROBLEMS
C
C DOUBLE PRECISION
C
C SUBROUTINE REQUIRED - INIVAL, CONT, DIFCOE, INTGND, DENERG, DCDP
C                       THESE ARE PROBLEM DEPENDENT
C                       WNFRQ, CHAREQ ARE ALSO PROBLEM DEPENDENT,
C                       BUT INCLUDED IN MAIN ROUTINES
C write(6,fmt) ---- write to screen
C
      REAL*8  RAD,TMASS,GS,RHO0,GA
      REAL*8  D,H,RHO,ELAMB,EMU,GR,XI,PHI,ETA
      COMMON/MODEL /NLAY,ISO,LAME,NDIV,DMAX,ILAY,
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000)
      REAL*8  RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,U,ENGY,DELTA
      REAL*8  A,F,Y
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6)
      REAL*8  YN,YB,SUM,Q
      COMMON/SOL   /YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1              WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2              AC(100)
C*    DIMENSION  FMT(20)
      CHARACTER*80 FMT
      DIMENSION  PERD(100),CMN(100),DC(100),CMX(100),DEP(100),EPSI(100),
     1           JYI(100),JDERI(100)
C
C character*40 infn
C	write(*,*) 'please input the input modle file name'
C	read(*,'(a)') infn
	open(5,file = 'scalifornia_refer.dat')
      READ(5,1)  ISET,JSET,KTIME                              !!!!!!!!!ISET,JSET,KTIME = 0
    1 FORMAT(3I5)
C ISET=INPUT DATA SET,  JSET=OUTPUT DATA SET
      IF(ISET.GT.7)  REWIND  ISET
c      WRITE(6,601)
      IF(JSET.GT.7)  REWIND  JSET
C*    CALL  INTIME(.TRUE.,.TRUE.,.FALSE.,.FALSE.,'  DISP  ',ITIME)
C
C
      READ(5,10)  MJOB
      
   10 FORMAT(I5)
C
C ************************************MJOB(outer loop)**************************************
C
      DO  9000  MJ=1,MJOB                                           !!!MJOB=1
c     WRITE(6,901)
  901 FORMAT(1H1)

      OPEN(UNIT=7,FILE='DERIV.DATA')
      OPEN(UNIT=8,FILE='Y.DATA')
C*    IF(KTIME.NE.0)
C*   1CALL  INTIME(.FALSE.,.FALSE.,.FALSE.,.TRUE.,'  DISP  ',ITIME)
C
      
      CALL RDMDL(ISET,JSET,IPMLS) 
      
                                                   !!!!RDMDL reads in earth model
      RMAX=RAD
C
c      WRITE(6,910)  NAME,NDIV,DMAX,RAD,NLAY,LAME,ISO,TMASS,GS,RHO0
  910 FORMAT(1H0,10H**********,20A4,10H**********/
     1       1H0,5HDIV.=,I5,2X,5HUP TO,F10.3/
     2       1H ,5HR   =,F9.3/
     3       1H ,5HLAY.=,I5/1H ,5HLAME=,I5/1H ,5HISO =,I5/
     4       1H ,5HMASS=,D15.5/1H ,5HGS  =,F11.5/1H ,5HRHO =,F11.5)
C
c      IPMLS=0        !!!! to print out model
      IF(IPMLS) 2341,2342,2341
2342  CALL MLIST
       
2341  READ(5,920) NJOB                                              !!!NJOB=1
      
  920 FORMAT(I5)
      NY=0
      ND=0
      IERROR=0
C
C***********************************NJOB(inner loop)******************************************
C
      DO  8000  NJ=1,NJOB
      READ(5,800) NP,DEPTHO,EPSO,MODE,JYO,JDERO,DDMNP,DDMXP
800   FORMAT(I5,F10.0,E10.3,3I2,F9.0,F10.0)
C NP=NUMBER OF POINTS IN ONE NJOB
C DEPTHO,EPSO,JYO,JDERO=DEFAULT VALUES FOR DEPTH,EPSIL,JY,JD
      READ(5,802)  FMT
C*802 FORMAT(20A4)
 802  FORMAT(A80)
      READ(5,FMT)  (PERD(NT),CMN(NT),DC(NT),CMX(NT),DEP(NT),EPSI(NT),
     1              JYI(NT),JDERI(NT),NT=1,NP)
	
      LY=0
      LD=0
      IF(IERROR.GT.1)  GO TO  8000
C
C PJOB
C
      DO  7000  NT=1,NP
      IF(IERROR.LE.1)  IERROR=0
      DEPTH=DEP(NT)
      EPSIL=EPSI(NT)
      JY   =JYI(NT)
      JD   =JDERI(NT)
C JY=PARAMETER FOR YLIST,  JD=PARAMETER FOR DLIST
      IF(DEPTH.EQ.0.0)  DEPTH=DEPTHO
      IF(EPSIL.EQ.0.0)  EPSIL=EPSO
      IF(JY.EQ.0)       JY=JYO
      IF(JD.EQ.0)       JD=JDERO
C
c      WRITE(6,700)
  700 FORMAT(1H1)
C*    IF(KTIME.EQ.0)  GO TO  701
C*    CALL  INTIME(.TRUE.,.TRUE.,.FALSE.,.FALSE.,'  PJOB  ',JTIME)
C*    CALL  INTIME(.FALSE.,.TRUE.,.FALSE.,.TRUE.,'  PJOB  ',JTIME)
  701 CONTINUE
       
       
      CALL  TTL1(MODE)
   
C
      CALL  TABLE(IBOTM,DEPTH,4)
      
      IF(IBOTM.LE.1)  GO TO  7000
      DEPTH=D(IBOTM)
C
      IF(PERD(NT))  730,730,710
  710 CONTINUE
      IF(MODE-5)  720,730,730
  720 CONTINUE
c      WRITE(6,721)  PERD(NT),DEPTH
  721 FORMAT(1H0,2HT=,F10.5,3X,6HDEPTH=,F8.3/1H0,5X,1HC,14X,5HDELTA)
      GO TO  750
  730 CONTINUE
      APERD=ABS(PERD(NT))
c      WRITE(6,731)  APERD,DEPTH
  731 FORMAT(1H0,2HN=,E12.5,3X,6HDEPTH=,F8.3/1H0,5X,1HT,14X,5HDELTA)
  750 CONTINUE
C
C TO FIND A ROOT OF CHARACTERISTIC EQUATION
C
      ITOP=1
      NMAX=0
       
      CALL INTERP(XX,PERD(NT),CMN(NT),DC(NT),CMX(NT),EPSIL,JUMP)
      
C*    IF(KTIME.NE.0)
C*   1CALL  INTIME(.FALSE.,.FALSE.,.TRUE.,.TRUE.,'  PJOB  ',JTIME)
C
      IF(JUMP)  6000,7100,6000
 7100 CONTINUE
      IF(IERROR-2)  7900,8500,9900
C
C TO COMPUTE EIGENFUNCTIONS
C
 6000 CONTINUE
      J=1
      IF((MODE.EQ.2).OR.(MODE.EQ.4))  CALL  INTEG
      J=2
      LY=LY+1
      DEP(LY)=DEPTH
C
      CALL  INTEG
      
C PRINT AND/OR PUNCH EIGENFUNCTIONS
C
      CALL YLIST(JY,JSET,DDMNP,DDMXP)
      
c      WRITE(6,601)
  601 FORMAT(1H0)
C*    IF(KTIME.NE.0)
C*   1CALL  INTIME(.FALSE.,.FALSE.,.TRUE.,.TRUE.,'  PJOB  ',JTIME)
C
C DERIVATIVES W.R.T. PARAMETERS
C
      IF(JD.EQ.0)  GO TO  7900
C
      CALL  DERIV
       
      LD=LD+1
C
C PRINT AND/OR PUNCH DERIVATIVES
      
      CALL DLIST(JD,JSET,DDMNP,DDMXP)
      
c      WRITE(6,601)
C*    IF(KTIME.NE.0)
C*   1CALL  INTIME(.FALSE.,.FALSE.,.TRUE.,.TRUE.,'  PJOB  ',JTIME)
C
 7900 CONTINUE
C
C END OF PJOB
C
 7000 CONTINUE
C
 8500 CONTINUE
      IF(LY.EQ.0)  GO TO  8900
      NY=NY+LY
      ND=ND+LD       
c      WRITE(6,901)                                                                  !!!!
c      WRITE(6,910)  NAME,NDIV,DMAX,RAD,NLAY,LAME,ISO,TMASS,GS,RHO0                  !!!!
      
      CALL  TTL1(MODE)
       
c      WRITE(6,851)                                                                 !!!!
  851 FORMAT(1H0,4X,4HK(N),13X,1HT,10X,1HC,10X,1HU,13X,4HK.E.,11X,
     1       5HY3/Y1,10X,5HERROR,8X,5HDEPTH)
c      WRITE(6,852)  (WNB(LYY),TT(LYY),CC(LYY),UU(LYY),ENG(LYY),ELL(LYY),            !!!!
c     1               AC(LYY),DEP(LYY),LYY=1,LY)                                     !!!!
	do lyy=1,ly
	write(9,*) tt(lyy),cc(lyy)
	enddo
  852 FORMAT(1H ,E14.7,F14.5,2F11.6,E16.7,2E15.5,F10.3)
C
 8900 CONTINUE
c      WRITE(6,601)
C*    IF(KTIME.NE.0)
C*   1CALL  INTIME(.FALSE.,.FALSE.,.TRUE.,.TRUE.,'  DISP  ',ITIME)
C
C ******************************END OF NJOB(inner loop)*******************************************
 8000 CONTINUE
C
 9100 CONTINUE
 
      IF(JSET.LE.7)  GO TO  9200
c      WRITE(6,930)  NAME
  930 FORMAT(1H0/1H0,20A4)
c      WRITE(6,931)  NY,ND
  931 FORMAT(1H ,I3,2X,24HSET(S) OF EIGENFUNCTIONS/
     1       1H ,I3,2X,21HSET(S) OF DERIVATIVES)
 9200 CONTINUE
C
C***************************** END OF MJOB(outer loop)*****************************************
C
 9000 CONTINUE
 9900 CONTINUE
      IF(ISET.GT.7)  REWIND  ISET
      IF(JSET.GT.7)  REWIND  JSET
c      WRITE(6,601)                                                       !!!!
C*    CALL  INTIME(.FALSE.,.FALSE.,.TRUE.,.TRUE.,'  DISP  ',ITIME)
      CLOSE(UNIT=7)
      CLOSE(UNIT=8)
	close(5)
	
	return
c      STOP
      END
