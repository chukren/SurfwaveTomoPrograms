      SUBROUTINE  DIFCOE
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
      REAL*8  RH,EL,EM,RVP,XC,PH,ET,R1,R2,EP
C
C COEFFICIENT MATRIX FOR DIFFERENTIAL EQUATIONS
C MANTLE WAVES
C ANISOTROPIC MODEL
C Y3+WNN*Y4=0
C
      DO  1  LL=1,6
      DO  1  LLL=1,6
      A(LL,LLL)=0.0
    1 CONTINUE
C
      RH=RHO(I)
      EL=ELAMB(I)
      EM=EMU(I)
      XC    =1.0D+0+XI(I)
      RVP   =EL+2.0D+0*EM
      PH    =1.0D+0+PHI(I)
      ET    =1.0D+0+ETA(I)
      EP=ET/PH
      R1    =1.0D+0/RC
      R2=R1*R1
C
      GO TO  (1000,1000,300,400,300,1000),MODE
C
C MANTLE LOVE WAVE
C TOROIDAL OSCILLATIONS
C
  300 CONTINUE
      IF(EM)  310,1000,310
  310 CONTINUE
      A(1,1)=2.0D+0*R1
      A(1,2)=1.0D+0/EM
      A(2,1)=-FRQ2*RH+R2*(WN2-2.0D+0)*XC*EM
      A(2,2)=-A(1,1)
      GO TO  1000
C
C MANTLE RAYLEIGH WAVES
C
  400 CONTINUE
      IF(EM)  410,450,410
  410 CONTINUE
      IF(J)  420,430,430
C
C DIFFERENTIAL EQUATION FOR CHARACTERISTIC EQUATION
C
  420 CONTINUE
      A(1,1)=R1*(3.0D+0-(2.0D+0*EP*EL)/RVP)
      A(1,4)=1.0D+0/EM
      A(1,5)=-1.0D+0/(PH*RVP)
      A(2,2)=-A(1,1)
      A(2,3)=4.0D+0*R2*(RVP-XC*EM-(ET*EP*EL*EL)/RVP)
      A(2,4)=-FRQ2*RH+A(2,3)
      A(2,5)=FRQ2*RH-R2*(WN2*(RVP-ET*EP*EL*EL/RVP)-2.0D+0*XC*EM)
      A(4,3)=-2.0D+0*R1*EP*EL/RVP
      A(4,4)=-R1*(1.0D+0+2.0D+0*EP*EL/RVP)
      A(5,3)=-2.0D+0*R1
      A(3,1)=-0.5D+0*WN2*A(2,3)
      A(3,4)=-0.5D+0*WN2*A(5,3)
      A(3,5)=-0.5D+0*WN2*A(4,3)
      A(4,1)=-A(2,5)
      A(4,2)=-A(1,5)
      A(5,1)=-A(2,4)
      A(5,2)=-A(1,4)
      A(5,5)=-A(4,4)
      GO TO  1000
C
C DIFFERENTIAL EQUATION FOR EIGENFUNCTIONS
C
  430 CONTINUE
      A(1,1)=R1*(1.0D+0-(2.0D+0*EP*EL)/RVP)
      A(1,2)=1.0D+0/(PH*RVP)
      A(3,1)=-R1
      A(3,3)=2.0D+0*R1
      A(3,4)=1.0D+0/EM
      A(4,1)=-2.0D+0*R2*(RVP-XC*EM-(ET*EP*EL*EL)/RVP)
      A(4,2)=-R1*EP*EL/RVP
      A(4,3)=-FRQ2*RH+R2*(WN2*(RVP-ET*EP*EL*EL/RVP)-2.0D+0*XC*EM)
      A(2,1)=-FRQ2*RH-2.0D+0*A(4,1)
      A(1,3)=-WN2*A(4,2)
      A(2,2)=-A(1,1)
      A(2,3)=WN2*A(4,1)
      A(2,4)=-WN2*A(3,1)
      A(4,4)=-A(3,3)
      GO TO  1000
C
C LIQUID LAYER
C
  450 CONTINUE
      A(1,1)=R1*(1.0D+0-2.0D+0*EP)
      A(1,2)=1.0D+0/(PH*EL)
      A(4,1)=-2.0D+0*R2*EL*(1.0D+0-ET*EP)
      A(4,2)=-R1*EP
      A(4,3)=-FRQ2*RH+R2*WN2*EL*(1.0D+0-ET*EP)
      A(1,3)=-WN2*A(4,2)
      A(2,1)=-FRQ2*RH-2.0D+0*A(4,1)
      A(2,3)=WN2*A(4,1)
C
      A(1,1)=A(1,1)-(A(1,3)*A(4,1))/A(4,3)
      A(1,2)=A(1,2)-(A(1,3)*A(4,2))/A(4,3)
      A(2,1)=A(2,1)-(A(2,3)*A(4,1))/A(4,3)
      A(2,2)=-A(1,1)
      GO TO  1000
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
