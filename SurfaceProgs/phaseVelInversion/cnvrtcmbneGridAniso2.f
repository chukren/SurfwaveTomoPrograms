c  prepare input for shear velocity inversion and aniso plotting
c  by combining aniso coefficients
c  from inversions at individual periods into single file with coefficients at 
c  all periods for each location.  
c  
c  combineGridAniso < combineGridAnisoInp
      parameter (npermax = 40, ngridmax = 20000)
      character*70 foutput,finput(npermax),ampout,circleout
      real*4 lon(ngridmax),lat(ngridmax),cs2th(npermax,ngridmax)
      real*4 stddvcs(npermax,ngridmax),sn2th(npermax,ngridmax)
      real*4 stddvsn(npermax,ngridmax)
      real*4 residcs(npermax,ngridmax),residsn(npermax,ngridmax)
      real*4 phip(3),maxdiff
      real*4 phi(npermax,ngridmax),amp1(npermax,ngridmax)
      real*4 stdphi(npermax,ngridmax),stdamp(npermax,ngridmax)
      integer function length
      pi = 3.14159
      conv2deg = 180./pi
c  foutput for combining all periods at same spot for later inversion for aniso(depth)
      read(*,*) foutput
      read(*,*) facmult, scalefac
      read(*,*) nper
      do i = 1, nper
        read(*,*) finput(i+10)
	open(i+10, file = finput(i+10))
	read(i+10,*) nxy
      enddo
      do i = 1,nper
        do ii = 1, nxy
          read(i+10,*) lon(ii),lat(ii),cs2th(i,ii),stddvcs(i,ii),
     1          residcs(i,ii), sn2th(i,ii), stddvsn(i,ii), residsn(i,ii)
          a1 = cs2th(i,ii)
	  a2 = sn2th(i,ii)
	  
	  stda1 = facmult*stddvcs(i,ii)
	  stda2 = facmult*stddvsn(i,ii)
c  a1 and a2 are cos2theta and sin2theta coefficients respectively
c  stda1 and stda2 are their standard deviations
c  covariance between a1 and a2 not considered here
c  phi is fast direction, theta azimuth from north
          phi(i,ii) = 0.5*atan2(a2,a1)
	  amp = sqrt(a1*a1+a2*a2)
	  amp1(i,ii) = amp
	  g11 = (-0.5*a2/amp**2)/stda1
	  g21 = (0.5*a1/amp**2)/stda2
	  g12 = a1/amp/stda1
c	  g12 = (0.5*a1/amp**2)/stda2
c	  g21 = a1/amp/stda1
	  g22 = a2/amp/stda2
c  construct components of gtransposeg where G is partial derivative matrix
	  gtga = g11*g11 + g21*g21
	  gtgb = g22*g22 + g12*g12
	  gtgc = g11*g12 + g21*g22
c  then taking advantage of analytical inverse of 2x2 symmetric matrix to get
c  diagonal elements that are variances of the model parameters, i.e., amp and
c  fast direction
c	  stdphi = sqrt(gtgb/(gtga*gtgb-gtgc*gtgc))
	  stdamp(i,ii) = sqrt(gtga/(gtga*gtgb-gtgc*gtgc))
c  non-linearity makes estimate of uncertainty in direction very unreliable
c  estimate crudely and conservatively by perturbing cos2theta and sin2theta
c  terms
c
c  perturb in opposite direction of sign of coefficient to make larger
c  angle perturbation
	  if (a1.ge.0.0) then
	    a1p = a1 - stda1
	    a1p2 = a1 - stda1/1.414
	  else
	    a1p = a1 + stda1
	    a1p2 = a1 + stda1/1.414
	  endif
	  if (a2.ge.0.0) then
	    a2p = a2 - stda2
	    a2p2 = a2 - stda2/1.414
	  else
	    a2p = a2 + stda2
	    a2p2 = a2 + stda2/1.414
	  endif
	  phip(1) = 0.5*atan2(a2,a1p)
	  phip(2) = 0.5*atan2(a2p,a1)
	  phip(3) = 0.5*atan2(a2p2,a1p2)
	  maxdiff = 0.0
	  do ip = 1,3
	    diff = abs(phip(ip)-phi(i,ii))
	    if (diff.gt.pi/2.0) diff = pi - diff
	    maxdiff = amax1(maxdiff,diff)
	  enddo
	  
	  phi(i,ii) = phi(i,ii)*conv2deg
	  stdphi(i,ii) = maxdiff*conv2deg
	  
	enddo
      enddo  
      open (50,file = foutput) 
      write(50,*) nxy
      write(50,*) nper   
      do ii = 1,nxy
        write(50,*) lon(ii),lat(ii)
        do i = 1,nper
	  write(50,*) cs2th(i,ii),stddvcs(i,ii),
     1          sn2th(i,ii), stddvsn(i,ii)
	enddo
      enddo
      do ii = 1,nxy
        write(50,*) lon(ii),lat(ii)
        do i = 1,nper
	  write(50,*) phi(i,ii),stdphi(i,ii),amp1(i,ii),
     1          stdamp(i,ii)
	enddo
      enddo
       
      do i = 1,nper
        close(i+10)
      enddo
      
      do i = 1, nper
	ampout = finput(i+10)
        lenth = length(ampout)
	ampout = ampout(1:lenth)//'amp'
	circleout = ampout(1:lenth)//'cir'
	open(2, file = circleout)
	open(3, file = ampout)
	do j = 1, nxy
	  diam = scalefac*stdamp(i,j)*20.0
	  slength = scalefac*amp1(i,j)*10.0
	  write(2,*) lon(j),lat(j),diam
	  write(3,*) lon(j),lat(j),phi(i,j),slength
	enddo
	close (2)
	close (3)
       enddo
      
      
      close(50)
      
      end
      
      function length (string)
      character string*(*)
      length = len(string)
      do while (string(length:length).eq.' ')
        length = length -1
      enddo
      return 
      end	    
   
