      SUBROUTINE  ENERGY
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
      REAL*8  SAVE,W1,W2
C
C PURPOSE -   TO CALCULATE ENERGY INTEGRALS USING THREE ORDINATES,
C     (I,N), (I-2,N-1), AND (I-4,N-2).
C             ON EXIT INPUT VALUES ARE RETURNED TO RC AND I.
C SUBROUTINE REQUIRED - INTGND
C
      SAVE=RC
      II=I
      H1=H(I)+H(I-1)
      H2=H(I-2)+H(I-3)
      W1=(H1+H2)/6.0D+0
      ALP=H2/H1
C
      W2    =W1*(2.0D+0-ALP)
      RC=RAD-D(I)
      DO  10  LL=1,L
      Y(LL)=YN(LL,N)
   10 CONTINUE
C
      CALL  INTGND
C
      DO  100  IS=1,ISUM
      SUM(IS)=SUM(IS)+W2*F(IS)
  100 CONTINUE
C
      W2    =W1*(2.0D+0+ALP+1.0D+0/ALP)
      I=I-2
      RC=RAD-D(I)
      DO  20  LL=1,L
      Y(LL)=YN(LL,N-1)
   20 CONTINUE
C
      CALL  INTGND
C
      DO  200  IS=1,ISUM
      SUM(IS)=SUM(IS)+W2*F(IS)
  200 CONTINUE
C
cc      W2    =W1*(2.0D+0-1.0D+0/ALP)
      W2    =W1*(-1.0D+0/ALP+2.0d0)
      I=I-2
      RC=RAD-D(I)
      DO  30  LL=1,L
      Y(LL)=YN(LL,N-2)
   30 CONTINUE
C
      CALL  INTGND
C
      DO  300  IS=1,ISUM
      SUM(IS)=SUM(IS)+W2*F(IS)
  300 CONTINUE
C
C EXIT
C
 1000 CONTINUE
      RC=SAVE
      I=II
      RETURN
      END
