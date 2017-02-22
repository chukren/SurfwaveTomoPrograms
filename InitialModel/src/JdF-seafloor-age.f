!======================================================================
! Read sea floor age from magnatic anomaly compiled by Wilson 
! Code is modified from Don Forsyth 
!
! Y.Y. Ruan  Brown  08/25/2014
!======================================================================
      subroutine crustage(xlon,ylat,cfloorage)
      real*4 polylat(50,50), polylon(50,50)
      real*4 maganomlat(50,20,50), maganomlon(50,20,50),magage(50,20)
      real*4 x(50),y(50), magint
      integer*4 npolypts(50), nmaganom(20), nmagpts(50,20)
      character*70 outmagfile
      logical in
      common /polymag/ polylat,polylon,maganomlat,maganomlon,magage,
     1    npoly, npolypts, nmaganom,nmagpts, outmagfile 

      outmagfile = 'PolygonsMagAnom'
      !read(*,*) outmagfile
      !read(*,*) xlon, ylat

      ! longitude needs to be -130 if west not +230
      ! read in polygons and magnetic anomalies
      call inputpolysmags

      ! find which polygon point lies within
      ipoly = 0
      in = .FALSE.

      do i = 1, npoly
        do j = 1, npolypts(i)
          x(j) = polylon(i,j)
          y(j) = polylat(i,j)
        enddo

        !call polygin(x,y,npolypts,xlon,ylat,in)
        call point_in_polygon(x,y,npolypts(i),xlon,ylat,in)
        if (in .eqv. .TRUE.) then
          ipoly = i
          go to 1010
        endif
      enddo

1010  continue

      if (in .eqv. .FALSE.) then
        ipoly  =  999
        magint = -999.0 
        cfloorage = magint
        go to 1020
      endif

      !write(*,*) "ipoly=",ipoly

      ! find nearest magnetic anomaly point within polygon
      dist1 = 1.0E+30
      do i = 1, nmaganom(ipoly)
        do j = 1, nmagpts(ipoly,i)
          dist = sqrt((maganomlon(ipoly,i,j)-xlon)**2
     &               +(maganomlat(ipoly,i,j)-ylat)**2)
          if (dist.le.dist1) then
            dist1 = dist
            magi1 = i
            magj1 = j
          endif
        enddo
      enddo

      !write(*,*) 'First point here:',magi1,magj1

      ! find second nearest magnetic anomaly point on same magnetic anomaly
      ! check only adjacent points since mag points digitized sequentially
      magi2 = magi1
      if (magj1.eq.1) then 
        magj2 = 2
      else 
        if (magj1.eq.nmagpts(ipoly,magi1)) then
          magj2 = magj1-1
        else
          dista = sqrt((maganomlon(ipoly,magi1,magj1-1)-xlon)**2
     &                +(maganomlat(ipoly,magi1,magj1-1)-ylat)**2)
          distb = sqrt((maganomlon(ipoly,magi1,magj1+1)-xlon)**2
     &                +(maganomlat(ipoly,magi1,magj1+1)-ylat)**2)
          if (dista.lt.distb) then
            magj2 = magj1-1
          else 
            magj2 = magj1+1
          endif
        endif
      endif

      !write(*,*) 'Second point here:',magi2,magj2

      ! find nearest magnetic anomaly point on a different magnetic anomaly
      dist2 = 1.0E+30
      do i = 1, nmaganom(ipoly)
        if (i.ne.magi1) then
        do j = 1, nmagpts(ipoly,i)
          dist = sqrt((maganomlon(ipoly,i,j)-xlon)**2
     &                +(maganomlat(ipoly,i,j)-ylat)**2)
          if (dist.le.dist2) then
            dist2 = dist
            magi3 = i
            magj3 = j
          endif
        enddo
        endif
      enddo

      !write(*,*) 'Third point here:',magi3,magj3

      !  given three magnetic anomaly points and their values, find plane through
      !  the points to be used for interpolation of age
      xa = maganomlon(ipoly,magi1,magj1)
      ya = maganomlat(ipoly,magi1,magj1)
      za = magage(ipoly,magi1)
      xb = maganomlon(ipoly,magi1,magj2)
      yb = maganomlat(ipoly,magi1,magj2)
      zb = magage(ipoly,magi1)
      xc = maganomlon(ipoly,magi3,magj3)
      yc = maganomlat(ipoly,magi3,magj3)
      zc = magage(ipoly,magi3)
      AB1 = xb-xa
      AB2 = yb-ya
      AB3 = zb-za
      AC1 = xc-xa
      AC2 = yc-ya
      AC3 = zc-za

      !===
      ! if cz ~ 0, point a, b, and c are almost on a same line,  interpolation 
      ! will blow off.
      !===
 
      ! see if the distance from point c to line ab is too close 
      dc = abs(AB2*xc - AB1*yc - xa*yb + xb*ya)
      dc = dc / sqrt(AB2*AB2 + AB1*AB1)

      if (dc .lt. 1.0E-01) then
        if (magj1.eq.1) then
          magj2 = 3
        else
          if (magj1.eq.nmagpts(ipoly,magi1)) then
            magj2 = magj1-2
          else
            dista = sqrt((maganomlon(ipoly,magi1,magj1-2)-xlon)**2
     &                 + (maganomlat(ipoly,magi1,magj1-2)-ylat)**2)
            distb = sqrt((maganomlon(ipoly,magi1,magj1+2)-xlon)**2
     &                 + (maganomlat(ipoly,magi1,magj1+2)-ylat)**2)
            if (dista.lt.distb) then
              magj2 = magj1-2
            else
              magj2 = magj1+2
            endif
          endif
        endif

        !write(*,*) 'Second point here:',magi2,magj2
        xb = maganomlon(ipoly,magi1,magj2)
        yb = maganomlat(ipoly,magi1,magj2)
        AB1 = xb-xa
        AB2 = yb-ya
        AB3 = zb-za
      endif

      cx =  AB2*AC3 - AC2*AB3
      cy = -AB1*AC3 + AB3*AC1
      cz =  AB1*AC2 - AB2*AC1
      const = -(cx*xa+cy*ya+cz*za)
      magint = -(cx*xlon + cy*ylat + const)/cz
  
      if (magint.lt.0.0) magint = 0.0     
      cfloorage = magint

      !write(*,*) ipoly, magage(ipoly,magi1),magage(ipoly,magi3)
1020  continue
      !write(*,*) magint
      
      return
      end
      
!=======================================================================
!  read polygon outlines of regions in which age is to be estimated
!  last point in each polygon should be same as first
!=======================================================================
      subroutine inputpolysmags
      real*4 polylat(50,50), polylon(50,50)
      real*4 maganomlat(50,20,50), maganomlon(50,20,50),magage(50,20)
      integer*4 npolypts(50), nmaganom(20), nmagpts(50,20)
      character*70 dummy, outmagfile

      common /polymag/ polylat,polylon,maganomlat,maganomlon,magage,
     &    npoly, npolypts, nmaganom,nmagpts, outmagfile 

      open(10, file = outmagfile)
      read(10,*) npoly
      do i = 1, npoly
        read(10,*) dummy
        read(10,*) polylon(i,1),polylat(i,1)
        do j = 2, 50
          read(10,*) polylon(i,j),polylat(i,j)
          !  polygon routine sometimes has a problem if successive points 
          !  have the same exact latitude - perturb slightly
          if (polylat(i,j).eq.polylat(i,j-1)) then
             polylat(i,j)=polylat(i,j)+.001
          endif
          if ((polylon(i,j).eq.polylon(i,1)).and.
     &        (polylat(i,j).eq.polylat(i,1))) then
            npolypts(i) = j
            go to 1000
          endif
        enddo
1000    continue
        read(10,*) nmaganom(i)
        do k = 1, nmaganom(i)
          read(10,*) magage(i,k),nmagpts(i,k)
          read(10,*)(maganomlon(i,k,m),maganomlat(i,k,m),m=1,nmagpts(i,k))
        enddo
      enddo
      close (unit= 10)
      return
      end

!=======================================================================
! test whether point (xtest,ytest) lies within a polygon
!  n (x,y) points must be supplied sequentially outlining polygon
!
!  in = 0 outside, in = 1, inside
!
!  find max and min y values of each line segment and of total set
!  this part could be spun off and supplied to routine so it would
!  only have to be done once, speeding up process
!=======================================================================
      subroutine polygin(x,y,n,xtest,ytest,in)
      real*4 x(50),y(50),xmax(50),ymax(50),xmin(50),ymin(50)
      in = 0
      xmaxx = -1.0E+30
      ymaxx = -1.0E+30
      xminn =  1.0E+30
      yminn =  1.0E+30

      do i = 1, n-1
        xmax(i) = amax1(x(i),x(i+1))
        ymax(i) = amax1(y(i),y(i+1))
        xmin(i) = amin1(x(i),x(i+1))
        ymin(i) = amin1(y(i),y(i+1))

        xmaxx = amax1(xmaxx,xmax(i))
        ymaxx = amax1(ymaxx,ymax(i))
        xminn = amin1(xminn,xmin(i))
        yminn = amin1(yminn,ymin(i))
      enddo

      xmax(n) = amax1(x(1),x(n))
      ymax(n) = amax1(y(1),y(n))
      xmin(n) = amin1(x(1),x(n))
      ymin(n) = amin1(y(1),y(n))

      if ((xtest.ge.xminn).and.(xtest.le.xmaxx)) then
        if ((ytest.ge.yminn).and.(ytest.le.ymaxx)) then
        !  shoot ray in positive y direction from (xtest,ytest).  If intersects
        !  even number (including zero) of bounding segments, point is outside
          inum = 0
	  do i = 1, n-1
	    if (x(i+1).ne.x(i)) then
	      yint = y(i)+(y(i+1)-y(i))*(xtest-x(i))/(x(i+1)-x(i))
	      if ((yint.ge.ymin(i)).and.(yint.le.ymax(i)).and.
     1             (yint.ge.ytest)) then
                inum=inum+1
	      endif
	    else
	      if ((xtest.eq.x(i)).and.((ytest.ge.ymin(i)).and.
     1               (ytest.le.ymax(i)))) then
     	        inum = inum+1
	      endif
	    endif
	  enddo

	  if (x(1).ne.x(n)) then
	    yint = y(n)+(y(1)-y(n))*(xtest-x(n))/(x(1)-x(n))
	    if ((yint.ge.ymin(n)).and.(yint.le.ymax(n)).and.
     1           (yint.ge.ytest)) then
              inum=inum+1
	    endif
	  else
	    if ((xtest.eq.x(n)).and.((ytest.ge.ymin(n)).and.
     1             (ytest.le.ymax(n)))) then
     	      inum = inum+1
	    endif
	  endif  

	  diff = float(inum)/2.0 - int(inum/2)
	  if (diff.gt.0.00001) in = 1
	endif
      endif
      return
      end



      subroutine point_in_polygon(x,y,n,xtest,ytest,isin)
      implicit none
      real, dimension(1:50) :: x, y
      real :: xtest, ytest, xint
      real :: xmax, xmin, ymax, ymin
      integer :: n, i, j
      logical :: isin

      isin = .FALSE.

      xmax = -1.0E+30
      ymax = -1.0E+30
      xmin =  1.0E+30
      ymin =  1.0E+30

      do i = 1, n
        xmax = max(xmax, x(i))
        xmin = min(xmin, x(i))
        ymax = max(ymax, y(i))
        ymin = min(ymin, y(i))
      enddo

      
      if ((xtest.ge.xmin).and.(xtest.le.xmax)) then
        if ((ytest.ge.ymin).and.(ytest.le.ymax)) then
          j = n-1
          do i = 1, n-1
            if((y(i).lt.ytest .and. y(j).ge.ytest) .or. (y(j).lt.ytest .and. y(i).ge.ytest)) then
              xint = x(i) + (ytest - y(i))/(y(j)-y(i))*(x(j)-x(i))
              if(xint .lt. xtest) then 
                isin = .not. isin 
              endif
            endif
            j = i
          enddo
        endif
      endif

      return
      end 
       
