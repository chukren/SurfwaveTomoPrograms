      SUBROUTINE  TTL1(MODE)
C     MODE = 4
      GO TO  (100,200,300,400,500,600),MODE
C
C      LOVE WAVE
       
  100 CONTINUE
      WRITE(6,101)
  101 FORMAT(1H0,22HLOVE WAVE (FLAT MODEL))
      GO TO  1000
C
C      RAYLEIGH WAVE
C
  200 CONTINUE
      
      WRITE(6,201)
  201 FORMAT(1H0,26HRAYLEIGH WAVE (FLAT MODEL))
      GO TO  1000
C
C      MANTLE LOVE WAVE
C
  300 CONTINUE
      WRITE(6,301)
  301 FORMAT(1H0,34HMANTLE LOVE WAVE (SPHERICAL MODEL))
      GO TO  1000
C
C      MANTLE RAYLEIGH WAVE
C
  400 CONTINUE
c      WRITE(6,401)
  401 FORMAT(1H0,38HMANTLE RAYLEIGH WAVE (SPHERICAL MODEL))
      
      GO TO  1000
      
C      TORSIONAL OSCILLATION
C
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(1H0,21HTORSIONAL OSCILLATION)
      GO TO  1000
C
C      SPHEROIDAL OSCILLATION
C
  600 CONTINUE
     
      WRITE(6,601)
  601 FORMAT(1H0,22HSPHEROIDAL OSCILLATION)
      
      GO TO  1000
      
 1000 CONTINUE
      
      RETURN
      END
