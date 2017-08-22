      SUBROUTINE  DCDP(DG)
      REAL*8  RAD,TMASS,GS,RHO0,GA                                    !!!!!RHO0  =DENSITY IN THE INNERMOST HOMOGENENEOUS CORE
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
      COMMON/SOL/  YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),         !!!!YN(1,*) YN(2,*) YN(3,*): density,dca(j,i),dcb(j,i)
     1             WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2             AC(100)
      REAL*8  DG,RH,EL,EM,RVP,XC,PH,ET,EP,DRO,DVP,DVS,DXC,DPH,DET,
     1        Y1,Y2,Y3,Y4,Y11,Y22,Y33,Y44,R2,CU,WW
C
C PARTIAL DERIVATIVES W.R.T.(with respect to) RHO, VP, VS, XI, PHI, AND ETA
C MANTLE WAVE (SPHERICAL MODEL)
C
      DRO=0.0
      DVP=0.0
      DVS=0.0
      DXC=0.0
      DPH=0.0
      DET=0.0
      RH=RHO(I)
      EL=ELAMB(I)
      EM=EMU(I)
      RVP   =EL+EM*2.0D+0
      XC    =1.0D+0+XI(I)
      PH    =1.0D+0+PHI(I)
      ET    =1.0D+0+ETA(I)
      EP=ET/PH
      Y1=YN(1,N)
      Y2=YN(2,N)
      Y3=YN(3,N)
      Y4=YN(4,N)
      Y11=Y1*Y1
      Y22=Y2*Y2
      R2=RC*RC
      CU    =0.5D+0*C/(ENGY*U)
C
      GO TO  (1000,1000,300,400,300,1000),MODE
C
C MANTLE LOVE WAVE (SPHERICAL MODEL)
C
  300 CONTINUE
      IF(EM.EQ.0.0)  GO TO  900
      DRO=-FRQ2*RH*Y11+Y22/EM+(WN2-2.0D+0)*XC*EM*Y11/R2
      DVS=2.0D+0*(Y22/EM+(WN2-2.0D+0)*XC*EM*Y11/R2)
      DXC=(WN2-2.0D+0)*XC*EM*Y11/R2
      GO TO  900
C
C MANTLE RAYLEIGH WAVE (SPHERICAL MODEL)
C
  400 CONTINUE
      IF(EM.lt.0.001)  GO TO  450
C        IF(EM.EQ.0.0)  GO TO  450
C
C SOLID LAYER
C
      WW    =2.0D+0*Y1-WN2*Y3
      Y33=Y3*Y3
      Y44=Y4*Y4
      DRO=-FRQ2*RH*(Y11+WN2*Y33)
     1   +(RVP-XC*EM-(EP*ET*EL*EL)/RVP)*WW*WW/R2
     2   +Y22/(PH*RVP)
     3   +WN2*(WN2-2.0D+0)*XC*EM*Y33/R2
     4   +WN2*Y44/EM
      DVP=2.0D+0*(((Y2+2.0D+0*ET*EM*WW/RC)**2)/(PH*RVP)
     1   +RVP*(1.0D+0-EP*ET)*WW*WW/R2)
      DVS=2.0D+0*(((-XC*EM+4.0D+0*EP*ET*EL*EM/RVP)*WW/R2
     1   -4.0D+0*EP*EM*Y2/(RVP*RC))*WW
     2   +WN2*(WN2-2.0D+0)*XC*EM*Y33/R2
     3   +WN2*Y44/EM)
      DXC=XC*EM*(-WW*WW+WN2*(WN2-2.0D+0)*Y33)/R2
      DPH=((Y2-ET*EL*WW/RC)**2)/(PH*RVP)
      DET=2.0D+0*EP*EL*WW*(Y2-ET*EL*WW/RC)/(RVP*RC)
      GO TO  900
C
C LIQUID LAYER
C
  450 CONTINUE
      A(4,3)=FRQ2*RH-WN2*EL*(1.0D+0-EP*ET)/R2
      A(4,1)=-2.0D+0*EL*(1.0D+0-EP*ET)/R2
      A(4,2)=-EP/RC
      Y3=(A(4,1)*Y1+A(4,2)*Y2)/A(4,3)
      Y33=Y3*Y3
      WW    =2.0D+0*Y1-WN2*Y3
      DRO=-FRQ2*RH*(Y11+WN2*Y33)
     1   +EL*(1.0D+0-EP*ET)*WW*WW/R2
     2   +Y22/(PH*EL)
      DVP=2.0D+0*(Y22/(PH*EL)+EL*(1.0D+0-EP*ET)*WW*WW/R2)
      DPH=((Y2-ET*EL*WW/RC)**2)/(PH*EL)
      DET=2.0D+0*EP*WW*(Y2-ET*EL*WW/RC)/RC
      
      GO TO  900
C
C
  900 CONTINUE
      YN(1,N)=CU*DRO
      YN(2,N)=CU*DVP
      YN(3,N)=CU*DVS
      YN(4,N)=CU*DXC
      YN(5,N)=CU*DPH
      YN(6,N)=CU*DET
      GO TO  1000

 1000 CONTINUE
      RETURN
      END
