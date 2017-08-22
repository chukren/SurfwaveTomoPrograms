      SUBROUTINE  GRAVTY
      REAL*8  RAD,TMASS,GS,RHO0,GA
      REAL*8  D,H,RHO,ELAMB,EMU,GR,XI,PHI,ETA
      COMMON/MODEL /NLAY,ISO,LAME,NDIV,DMAX,ILAY,
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000)
      REAL*8  F1,F2,F3,H1,H2,C1,C2,C3,EM,A
      REAL*8 PEI,G
C*    REAL*8  PEI/12.56637061435917/
      DATA PEI/12.56637061435917/
C*    REAL  G/6.67E-8/
      DATA G/6.67D-8/
C
C PURPOSE - TO COMPUTE THE GRAVITY DISTRIBUTION
C INPUT - TMASS=TOTAL MASS OF THE EARTH IN C.G.S.
C         GS   =SURFACE VALUE OF GRAVITY IN C.G.S.
C         RHO0 =DENSITY OF THE INNERMOST HOMOGENEOUS SPHERE
C IF(TMASS.EQ.0)  TMASS MAY BE CALCULATED FROM GS
C IF(GS   .EQ.0)  INNER SPHERE R.LT.R(NLAY) IS ASSUMED TO BE HOMOGENEOUS
C                 WITH DENSITY OF RHO0
C IF(RHO0 .EQ.0)  RHO(NLAY) IS ASSIGNED TO RHO0
C OUTPUT - GR  =GRAVITY IN (10*5) C.G.S.
C
      GA=PEI*G
      I=NLAY
      EM=0.0
      GR(I)=0.0
      F3=RHO(I)*((RAD-D(I))**2)
  100 CONTINUE
      IF(H(I))  200,110,200
  110 CONTINUE
      I=I-1
      GR(I)=EM
      F3=RHO(I)*((RAD-D(I))**2)
      IF(I-1)  500,500,200
  200 CONTINUE
      F1=F3
      H1=H(I)
      H2=H(I-1)
      A=H2/H1
      F2=RHO(I-1)*((RAD-D(I-1))**2)
      F3=RHO(I-2)*((RAD-D(I-2))**2)
      C1    =2.0D+0+A/(1.0D+0+A)
      C2    =3.0D+0+1.0D+0/A
      C3    =-1.0D+0/(A*(1.0D+0+A))
      GR(I-1)=GR(I)+(H1/6.0D+0)*(C1*F1+C2*F2+C3*F3)
      C1    =2.0D+0-A
      C2    =2.0D+0+A+1.0D+0/A
      C3    =2.0D+0-1.0D+0/A
      EM    =GR(I)+((H1+H2)/6.0D+0)*(C1*F1+C2*F2+C3*F3)
      GR(I-2)=EM
      I=I-2
      IF(I-1)  500,500,100
C
  500 CONTINUE
      IF(TMASS)  510,520,510
  510 CONTINUE
      EM    =(TMASS*1.0D-15)/PEI-GR(1)
      RHO0  =(3.0D+0*EM)/((RAD-D(NLAY))**3)
      GO TO  570
  520 CONTINUE
      IF(GS)  530,540,530
  530 CONTINUE
      TMASS =((PEI*GS*RAD*RAD)/GA)*1.0D+10
      GO TO  510
  540 CONTINUE
      IF(RHO0)  550,560,550
  550 CONTINUE
      EM=RHO0*((RAD-D(NLAY))**3)/3.0D+0
      TMASS =PEI*(EM+GR(1))*1.0D+15
      GO TO  570
  560 CONTINUE
      RHO0=RHO(NLAY)
      GO TO  550
  570 CONTINUE
      DO  600  I=1,NLAY
      GR(I)=(GA*(GR(I)+EM))/((RAD-D(I))**2)
  600 CONTINUE
      GS    =GR(1)*1.0D+5
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END