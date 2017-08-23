      SUBROUTINE  DERIV
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
      REAL*8  W1,W2,W3,DG1,DG2,DG3,F1,F2,F3
      REAL*8  DGINTG
C
C PURPOSE - TO CALCULATE DERIVATIVES OF EIGENVALUE
C SUBROUTINE REQUIRED - DGINTG, DCDP            !!!!!!!!!!!!!DCDP: PARTIAL DERIVATIVES W.R.T.(with respect to) RHO, VP, VS, XI, PHI, AND ETA
C
      DG3=0.0
      N=0
      I=ITOP-1
C
  100 CONTINUE
      N=N+1
      I=I+1
      RC=RAD-D(I)
      IF(MODE.EQ.6)  F3=DGINTG(RC,I,N)
      
      CALL  DCDP(DG3)
      
      IF(N.GE.NMAX)  GO TO  1000
      
  200 CONTINUE
	 
      IF(MODE.NE.6)  GO TO  250
C
      H1=H(I+4)+H(I+3)
      H2=H(I+2)+H(I+1)
      F2=DGINTG(RC-H2,I+2,N+1)
      F1=DGINTG(RC-H1-H2,I+4,N+2)
      ALP=H2/H1
      W1    =-(ALP*ALP)/(1.0D+0+ALP)
      W2    =3.0D+0+ALP
      W3    =2.0D+0+1.0D+0/(1.0D+0+ALP)
      DG2   =DG3+(H2/6.0D+0)*(W1*F1+W2*F2+W3*F3)
      W1    =2.0D+0-ALP
      W2    =2.0D+0+ALP+1.0D+0/ALP
      W3    =2.0D+0-1.0D+0/ALP
      DG1   =DG3+((H1+H2)/6.0D+0)*(W1*F1+W2*F2+W3*F3)
C
  250 CONTINUE
      I=I+2
      N=N+1
      RC=RAD-D(I)
 	
      CALL  DCDP(DG2)
      I=I+2
      N=N+1
      RC=RAD-D(I)

      CALL  DCDP(DG1)
C
  300 CONTINUE
      F3=F1
      DG3=DG1
      IF(N.GE.NMAX)  GO TO  1000
      IF(H(I+1))  200,100,200
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
C ---------------------------------------------------------
      REAL*8   FUNCTION  DGINTG(RR,II,NN)
C*    REAL     FUNCTION  DGINTG*8(RR,II,NN)
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
      REAL*8  RR,R2,RH,EL,PH,ET,EP,Y1,Y3,GC
C
C PURPOSE - TO CALCULATE THE INTEGRAND OF DG/DRHO
C
      RH=RHO(II)
      EL=ELAMB(II)
      PH    =1.0D+0+PHI(II)
      ET    =1.0D+0+ETA(II)
      EP=ET/PH
      GC=GR(II)
      R2=RR*RR
      Y1=YN(1,NN)
      Y3=YN(3,NN)
      
      IF(EMU(II).gt.0.001)  GO TO  100
C LIQUID LAYER
      Y3    =((RH*GC/RR-2.0D+0*EL*(1.0D+0-ET*EP)/R2)*Y1
     1      -EP*YN(2,NN)/RR-RH*Y3/RR)
     2      /(FRQ2*RH-WN2*EL*(1.0D+0-ET*EP)/R2)
C
  100 CONTINUE
      DGINTG=(GA*RH*(2.0D+0*Y1-WN2*Y3)*Y1)/(RR*R2)
C
      RETURN
      END
