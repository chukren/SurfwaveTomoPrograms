c  cnvrtaniso2.f  
c
c  just like cnvrtaniso.f except for node coefficients instead of regions
c
c  converts cos2theta and sin2theta terms into amplitude and direction of
c  splitting, plus calculates uncertainty based on standard deviations of
c  coefficients
      parameter (nmaxreg = 1000)
      character*35 outputfile,setdescription(100),regiondescp(nmaxreg)
      character*35 circleout
      real*4 phip(3),maxdiff
      real*4 phi(nmaxreg,nmaxreg),amp1(nmaxreg,nmaxreg)
      real*4 stdphi(nmaxreg,nmaxreg),stdamp(nmaxreg,nmaxreg)
      real*4 rlon(nmaxreg),rlat(nmaxreg)
      integer function length
      pi = 3.14159
      conv2deg = 180./pi
      read(*,*) outputfile
      open(10, file = outputfile)
      read(*,*) facmult, scalefac
c  Like for phase velocities, the nominal standard deviations are too small by
c  factor of two to three, based on period to period oscillations.  facmult
c  multiplies standard deviations by a factor of facmult to correct
c
c  scalefac is scale in inches/0.1km/s for aniso vector on GMT plot
      read(*,*) nregions
      read(*,*) nsets
      write(10,*) nsets
      write(10,*) 'phi    stdphi   amp     stadamp'
      do i = 1, nregions
c  approximate coordinates for plotting aniso vectors on map
        read(*,*) rlon(i),rlat(i),kk, regiondescp(i)
c	write(*,*) i, rlon(i),rlat(i)
      enddo
      do j = 1, nsets
        read(*,*) numpairs, setdescription(j)
	write(10,*) numpairs,'  ', setdescription(j)
        do i = 1, numpairs
	  read(*,*) regiondescp(i)
c	  write(*,*) setdescription(j),regiondescp(i)
	  read(*,*) dummy
	  read(*,*) a1,stda1,a2,stda2
	  stda1 = facmult*stda1
	  stda2 = facmult*stda2
c  a1 and a2 are cos2theta and sin2theta coefficients respectively
c  stda1 and stda2 are their standard deviations
c  covariance between a1 and a2 not considered here
c  phi is fast direction, theta azimuth from north
          phi(i,j) = 0.5*atan2(a2,a1)
	  amp = sqrt(a1*a1+a2*a2)
	  amp1(i,j) = amp
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
	  stdamp(i,j) = sqrt(gtga/(gtga*gtgb-gtgc*gtgc))
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
	    diff = abs(phip(ip)-phi(i,j))
	    if (diff.gt.pi/2.0) diff = pi - diff
	    maxdiff = amax1(maxdiff,diff)
	  enddo
	  
	  phi(i,j) = phi(i,j)*conv2deg
	  stdphi(i,j) = maxdiff*conv2deg
	  
	  write(10,*) regiondescp(i)
	  write(10,*) phi(i,j),stdphi(i,j),amp1(i,j),stdamp(i,j)
	enddo
      enddo
      do i = 1, numpairs
        write(10,*) regiondescp(i)
        do j = 1, nsets
           write(10,*) setdescription(j),
     1     phi(i,j),stdphi(i,j),amp1(i,j),stdamp(i,j)
        enddo
      enddo
c  set up files for plotting aniso on phase velocity maps
c  file 11 for vectors.  Indicate uncertainties by drawing circles
c  at ends of vectors with radius = stdamp
      do j = 1, nsets
        open (11, file = setdescription(j))
	do i = 1, numpairs
	  slength = scalefac*amp1(i,j)*10.0
	  write(11,*) rlon(i),rlat(i),phi(i,j),slength
	enddo
	close (11)
      enddo
      do j = 1, nsets
	lenth = length(setdescription(j))
	circleout =setdescription(j)
	circleout = circleout(1:lenth)//'c'
	open(12, file = circleout)
	do i = 1, numpairs
	  diam = scalefac*stdamp(i,j)*20.0
	  write(12,*) rlon(i),rlat(i),diam
	enddo
	close (12)
      enddo
 	
      close (10)
      end
      
      function length (string)
      character string*(*)
      length = len(string)
      do while (string(length:length).eq.' ')
        length = length -1
      enddo
      return 
      end	    
          
