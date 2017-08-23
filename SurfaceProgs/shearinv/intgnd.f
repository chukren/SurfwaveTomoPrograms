      SUBROUTINE  INTGND
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
      REAL*8  RH,EL,EM,RVP,XC,PH,ET,EP,R2,Y3,Y33,Y44,W,WW
C
C TO COMPUTE INTEGRANDS OF ENERGY INTEGRALS
C MANTLE WAVE (SPHERICAL MODEL)
C
      DO  1  IS=1,ISUM
      F(IS)=0.0
    1 CONTINUE
C
      RH=RHO(I)
      EL=ELAMB(I)
      EM=EMU(I)
      XC    =1.0D+0+XI(I)
      PH    =1.0D+0+PHI(I)
      ET    =1.0D+0+ETA(I)
      RVP   =EL+EM*2.0D+0
      EP=ET/PH
      R2=RC*RC
C
      GO TO  (1000,1000,300,400,300,1000),MODE
C
C MANTLE LOVE WAVE (SPHERICAL MODEL)
C
  300 CONTINUE
      IF(EM.EQ.0.0)  GO TO  1000
      F(1)=RH*Y(1)*Y(1)
      F(3)=XC*EM*((Y(1)/RC)**2)
      F(2)=(WN2-2.0D+0)*F(3)+(Y(2)*Y(2))/EM
      F(3)=RAD*RAD*F(3)
      GO TO  1000
C
C MANTLE RAYLEIGH WAVE (SPHERICAL MODEL)
C
  400 CONTINUE
      IF(EM.EQ.0.0)  GO TO  450
C
C SOLID LAYER
C
      W=RVP-EP*ET*EL*EL/RVP
      WW    =2.0D+0*Y(1)-WN2*Y(3)
      Y33=Y(3)*Y(3)
      Y44=Y(4)*Y(4)
      F(1)=RH*(Y(1)*Y(1)+WN2*Y33)
      F(2)=(W-XC*EM)*WW*WW/R2
     1    +Y(2)*Y(2)/(PH*RVP)
     2    +WN2*(WN2-2.0D+0)*XC*EM*Y33/R2
     3    +WN2*Y44/EM
      F(3)=-FRQ2*RH*Y33
     1    -2.0D+0*W*Y(3)*WW/R2
     2    -2.0D+0*EP*EL*Y(2)*Y(3)/(RVP*RC)
     3    +2.0D+0*XC*EM*Y(3)*(2.0D+0*Y(1)-Y(3))/R2
     4    +Y44/EM
      F(3)=RAD*RAD*F(3)
      GO TO  1000
C
C LIQUID LAYER
C
  450 CONTINUE
      W     =EL*(1.0D+0-EP*ET)
      A(4,3)=FRQ2*RH-WN2*W/R2
      A(4,1)=-2.0D+0*W/R2
      A(4,2)=-EP/RC
      Y3=(A(4,1)*Y(1)+A(4,2)*Y(2))/A(4,3)
      WW    =2.0D+0*Y(1)-WN2*Y3
      Y33=Y3*Y3
      F(1)=RH*(Y(1)*Y(1)+WN2*Y33)
      F(2)=W*WW*WW/R2
     1    +Y(2)*Y(2)/(PH*EL)
      F(3)=-FRQ2*RH*Y33
     1    -2.0D+0*W*Y3*WW/R2
     2    -2.0D+0*EP*Y(2)*Y3/RC
      F(3)=RAD*RAD*F(3)
      GO TO  1000
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
