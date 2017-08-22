      SUBROUTINE  SPHERF(PHI1,PHI2,PSI1,PSI2,FN,X,M)
      REAL*8  PHI1,PHI2,PSI1,PSI2,FN,FNN,X,A,AA,B,BB,C
C
      PHI1=0.0
      PHI2=0.0
      PSI1=0.0
      PSI2=0.0
      A=1.0D+0
      B=1.0D+0
      AA=1.0D+0
      BB=1.0D+0
      C=2.0D+0
      FNN   =2.0D+0*FN
C
      DO  100  MM=1,M
      PHI1=PHI1+A
      PHI2=PHI2+AA
      PSI1=PSI1+B
      PSI2=PSI2+BB
      E=(A*AA*B*BB)/(PHI1*PHI2*PSI1*PSI2)
      IF(ABS(E).LE.1.0E-20)  GO TO  200
      A     =-(A*X)/(C*(FNN+C+1.0D+0))
      AA   =-(AA*X)/(C*(FNN+C+3.0D+0))
      B     =-(B*X)/((C+2.0D+0)*(FNN+C+3.0D+0))
      BB    =-(BB*X)/((C+2.0D+0)*(FNN+C+5.0D+0))
      C=C+2.0D+0
  100 CONTINUE
c      WRITE(6,101)  X,FN,M
  101 FORMAT(1H0,11HARGUMENT X=D12.5,3X,12HIS TOO LARGE,3X,2HN=F4.0,3X,
     1       2HM=I3,3X,8H(SPHERF)//)
  200 CONTINUE
      RETURN
      END
