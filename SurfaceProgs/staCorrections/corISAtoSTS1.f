c corISAtoSTS1.f
c  use in conjuntion with sac macro brsortcut.m to correct misbehaving ISA 
c  to equivalent of STS1 vertical response

c  Empirical correction - no significance to functional form of correction
c  which are just convenient, smooth fits to observed spectral transfer
c  functions.  CAUTION  - this should be believed only between period of 
c  16 s and 143 s where we have inversions for station corrections based
c  on uniform phase velocity and wave parameters


      character*80 f1a, f1p, f2a, f2p
      parameter (nmax = 655360)
      dimension x1(nmax), x2(nmax), x3(nmax), x4(nmax)      
      f1a = 'temp.am'
      f2a = 'tempf.am'
      f1p = 'temp.ph'
      f2p = 'tempf.ph' 
      pi = 3.1415928
      twopi = 2.0*pi
      tenlg = 2.302585
      
      call rsac1(f1a,x1,npts,beg,delf,nmax,nerr)
      call rsac1(f1p,x2,npts,beg,delf,nmax,nerr)
      if (npts.gt.nmax) then
         write(*,*)  'time series too long'
      end if
      do i = 2, npts       
        freq= delf*(i-1)
	flog = alog10(freq)
	x3(i) = x1(i)/(0.55*flog +2.014)/0.75
	if (freq.gt.0.065) x3(i) = x1(i)/1.3/0.75
        x4(i) = x2(i) - 0.0166*twopi
	if (freq.le. 0.0113) then
 	  x4(i) = x2(i) + (-flog*0.0387 - 0.1081)*twopi
        endif
	if ((freq.gt.0.013).and.(freq.le.0.0231)) then
	  x4(i) = x2(i) + (flog*0.0521 + 0.0686)*twopi
	endif
      enddo
      x3(1) = 0.0
      x4(1) = 0.0
	call wsac0(f2a,xdummy,x3,nerr)
	call wsac0(f2p,xdummy,x4,nerr)
        close(unit=10)
	stop	
      end      
