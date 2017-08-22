      SUBROUTINE  DIVIDE
      REAL*8  RAD,TMASS,GS,RHO0,GA
      REAL*8  D,H,RHO,ELAMB,EMU,GR,XI,PHI,ETA
      COMMON/MODEL /NLAY,ISO,LAME,NDIV,DMAX,ILAY,
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000)
      REAL*8  D1,D2,D3,Q1,Q2,Q3,Q4,Q5,Q6,AR,BR,CR,AL,BL,CL,AM,BM,CM,
     1        AX,BX,CX,AP,BP,CP,AE,BE,CE,X,DX,DD,HH
C
C PURPOSE - TO DIVIDE ONE STEP INTO NDIV STEPS
C           LAYERS BELOW DMAX REMAIN INTACT
C
       
      N=NDIV-1
      DD=0.0
      I=0
      J=0
      IF(NDIV-1)  500,500,10
   10 CONTINUE
      DD=HI(1)
      HH=DD
C*    DX    =1.0D+0/DFLOAT(NDIV)
      DX    =1.0D+0/DBLE(NDIV)
C
  100 CONTINUE
      I=I+1
      J=J+1
      D(J)=DD
      H(J)=HH
      RHO(J)=RHOI(I)
      ELAMB(J)=VPI(I)
      EMU(J)=VSI(I)
      XI(J)=XII(I)
      PHI(J)=PHII(I)
      ETA(J)=ETAI(I)
      IF(I-ILAY)  110,1000,1000
  110 CONTINUE
   
      IF(DD-DMAX)  120,500,500
C
C BACKWARD INTERPOLATION
C
  120 CONTINUE
      D1=HI(I+1)
      IF(I+1-ILAY)  130,400,400
  130 CONTINUE
      D2=HI(I+2)
      IF(D2)  140,400,140
  140 CONTINUE
C
      Q1=D1/(D1+D2)
      Q2=-D1/D2
      Q3=-Q1*Q2
      Q4    =-1.0D+0-Q1
      Q5    =1.0D+0-Q2
      Q6=-Q3
C
      AR=Q1*RHOI(I)+Q2*RHOI(I+1)+Q3*RHOI(I+2)
      BR=Q4*RHOI(I)+Q5*RHOI(I+1)+Q6*RHOI(I+2)
      CR=RHOI(I)
      AL=Q1*VPI(I)+Q2*VPI(I+1)+Q3*VPI(I+2)
      BL=Q4*VPI(I)+Q5*VPI(I+1)+Q6*VPI(I+2)
      CL=VPI(I)
      AM=Q1*VSI(I)+Q2*VSI(I+1)+Q3*VSI(I+2)
      BM=Q4*VSI(I)+Q5*VSI(I+1)+Q6*VSI(I+2)
      CM=VSI(I)
      AX=Q1*XII(I)+Q2*XII(I+1)+Q3*XII(I+2)
      BX=Q4*XII(I)+Q5*XII(I+1)+Q6*XII(I+2)
      CX=XII(I)
      AP=Q1*PHII(I)+Q2*PHII(I+1)+Q3*PHII(I+2)
      BP=Q4*PHII(I)+Q5*PHII(I+1)+Q6*PHII(I+2)
      CP=PHII(I)
      AE=Q1*ETAI(I)+Q2*ETAI(I+1)+Q3*ETAI(I+2)
      BE=Q4*ETAI(I)+Q5*ETAI(I+1)+Q6*ETAI(I+2)
      CE=ETAI(I)
      D3=D1
C
  200 CONTINUE
      X=DX
      HH=D3*DX
C
      DO  300  NN=1,N
      J=J+1
      DD=DD+HH
      D(J)=DD
      H(J)=HH
      RHO(J)  =X*(AR*X+BR)+CR
c      write(*,111)RHO(J),ELAMB(J),EMU(J),XI(J),PHI(J),ETA(J) 
c111	format(6F10.2)       
      ELAMB(J)=X*(AL*X+BL)+CL
      EMU(J)  =X*(AM*X+BM)+CM
      XI(J)   =X*(AX*X+BX)+CX
      PHI(J)  =X*(AP*X+BP)+CP
      ETA(J)  =X*(AE*X+BE)+CE
      X=X+DX
  300 CONTINUE
      
      I=I+1
      J=J+1
      DD=DD+HH
      D(J)=DD
      H(J)=HH
      RHO(J)=RHOI(I)
      ELAMB(J)=VPI(I)
      EMU(J)=VSI(I)
      XI(J)=XII(I)
      PHI(J)=PHII(I)
      ETA(J)=ETAI(I)
C
      IF(I-ILAY)  310,1000,1000
  310 CONTINUE
      IF(DD-DMAX)  320,500,500
  320 CONTINUE
      HH=HI(I+1)
      IF(HH)  330,100,330
  330 CONTINUE
C
C FORWARD INTERPOLATION
C
      D1=HI(I)
      IF(D1)  340,400,340
  340 CONTINUE
      D2=HI(I+1)
C
      Q2=-D2/D1
      Q3=D2/(D1+D2)
      Q1=-Q2*Q3
      Q4=-Q1
      Q5    =-1.0D+0-Q2
      Q6    =1.0D+0-Q3
C
      AR=Q1*RHOI(I-1)+Q2*RHOI(I)+Q3*RHOI(I+1)
      BR=Q4*RHOI(I-1)+Q5*RHOI(I)+Q6*RHOI(I+1)
      CR=RHOI(I)
      AL=Q1*VPI(I-1)+Q2*VPI(I)+Q3*VPI(I+1)
      BL=Q4*VPI(I-1)+Q5*VPI(I)+Q6*VPI(I+1)
      CL=VPI(I)
      AM=Q1*VSI(I-1)+Q2*VSI(I)+Q3*VSI(I+1)
      BM=Q4*VSI(I-1)+Q5*VSI(I)+Q6*VSI(I+1)
      CM=VSI(I)
      AX=Q1*XII(I-1)+Q2*XII(I)+Q3*XII(I+1)
      BX=Q4*XII(I-1)+Q5*XII(I)+Q6*XII(I+1)
      CX=XII(I)
      AP=Q1*PHII(I-1)+Q2*PHII(I)+Q3*PHII(I+1)
      BP=Q4*PHII(I-1)+Q5*PHII(I)+Q6*PHII(I+1)
      CP=PHII(I)
      AE=Q1*ETAI(I-1)+Q2*ETAI(I)+Q3*ETAI(I+1)
      BE=Q4*ETAI(I-1)+Q5*ETAI(I)+Q6*ETAI(I+1)
      CE=ETAI(I)
      D3=D2
      GO TO  200
C
C LINEAR INTERPOLATION
C
  400 CONTINUE
      AR=0.0
      BR=RHOI(I+1)-RHOI(I)
      CR=RHOI(I)
      AL=0.0
      BL=VPI(I+1)-VPI(I)
      CL=VPI(I)
      AM=0.0
      BM=VSI(I+1)-VSI(I)
      CM=VSI(I)
      AX=0.0
      BX=XII(I+1)-XII(I)
      CX=XII(I)
      AP=0.0
      BP=PHII(I+1)-PHII(I)
      CP=PHII(I)
      AE=0.0
      BE=ETAI(I+1)-ETAI(I)
      CE=ETAI(I)
      D3=HI(I+1)
      GO TO  200
C
  500 CONTINUE
      I=I+1
      DO  600  II=I,ILAY
      J=J+1
      HH=HI(II)
      DD=DD+HH
      D(J)=DD
      H(J)=HH
      RHO(J)=RHOI(II)
      ELAMB(J)=VPI(II)
      EMU(J)=VSI(II)
      XI(J)=XII(II)
      PHI(J)=PHII(II)
      ETA(J)=ETAI(II)
  600 CONTINUE
C
C EXIT
C
 1000 CONTINUE
      NLAY=J
      RETURN
      END
