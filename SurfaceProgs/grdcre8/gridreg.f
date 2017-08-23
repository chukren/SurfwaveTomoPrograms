***************************************************************************
*   Create regularly spaced grid points in an area with center point (slat,slon). *
*   In order to avoid distortion at high latitude, change coordinates so   *
*   center of area lies along the equator of coordinate system and there is 
*   compensation for curvature.                       *
*************************************************************************** 
c  gridreg < brwestreginp.dat
       parameter (nxptmx=100, nyptmx=100)
       real*4 dazim(nxptmx),ddelt(nyptmx),glat(nxptmx,nyptmx)
       real*4 glon(nxptmx,nyptmx)
       character*70 fdelta,fazim,foutput
       radius = 6371.
       pi = 3.14158
       convdeg = 2.0*pi/360.
c   'Enter the center point (slat, slon):'
       read(*,*) slat,slon
       sdelta=90.
c      sazimz is azimuth from the center point that will tilt grid relative to North
       read(*,*) sazimz
c  delx is increment in degrees at equator - the true increment will vary with latitude from
c  equator of projection
       read(*,*) nxpt, delx,begx
c to follow our traditional convention of listing points from right to left and bottom to top.
c  the degree increment for distance from the pole, dely, should be negative and the
c  beginning point should be farthest in degrees from projection pole (assuming sazimz is
c  northwards), (begy = 6 (degrees south of center point) is 96 degrees from pole)
c  begx would be negative and increase with positive delx
       read(*,*) nypt,dely, begy
       do i = 1, 4
         read(*,*) dummy, dummy2
c  dummy points are the box outline points to be used in inversion program
       enddo
       read(*,*) foutput
	open(3, file=foutput, status='unknown')
       degdist = 2.0*pi*radius/360.
       dx = delx*degdist
       dy = dely*degdist
       call gohead(slat,slon,sdelta,sazimz,plat,plon)
       call disthead(plat,plon,slat,slon,delta,azimz)
       ntot = nxpt*nypt
       write(3,*) ntot,nxpt, nypt,dx,dy
       do 150 j=1,nypt
         delyinc = (j-1)*dely +begy	 
         delt = delta + delyinc	 
         do 100 i=1,nxpt
	   azim = azimz + (begx + (i-1)*delx)/cos(convdeg*delyinc)
	   call gohead(plat,plon,delt,azim,glat(i,j),glon(i,j))
	   if (glon(i,j).gt.180.0) then
	     glon(i,j) = glon(i,j) - 360.0
	   endif
  100	 enddo
  150  enddo 
       do i = 1, nxpt
         do j = 1, nypt
	   write(3,200) glat(i,j),glon(i,j)
	 enddo
       enddo
       close(3)  
  200  format(f7.3,f10.3)
       end  
       
       
       
      subroutine disthead(slat,slon,flat,flon,delta,azim)
c  Calculates distance and azimuth on sphere from starting point s 
c  to finishing point f
      dtor= 3.1415928/180.
      slt = slat*dtor
      sln = slon*dtor
      flt = flat*dtor
      fln = flon*dtor
      delta = acos(sin(slt)*sin(flt)+cos(slt)*cos(flt)*cos(fln-sln))
      azim = atan2(sin(fln-sln)*cos(flt),
     *  sin(flt)*cos(slt) - cos(fln-sln)*cos(flt)*sin(slt))
      delta = delta/dtor
      azim = azim/dtor
      return
      end

      subroutine gohead(slat,slon,delta,azim,flat,flon)
c Calculates final latitude and longitude f when starting at point s
c traveling delta degrees on sphere at initial heading azim
  
      dtor= 3.1415928/180.
      slt = slat*dtor
      dlta = delta*dtor
      azm = azim*dtor
      flat = asin(cos(slt)*sin(dlta)*cos(azm) + sin(slt)*cos(dlta))
      flon = atan2(sin(dlta)*sin(azm), 
     *   cos(slt)*cos(dlta) - sin(slt)*sin(dlta)*cos(azm))
      flat = flat/dtor
      flon = slon + flon/dtor
      if (flon.gt.360.) flon = flon - 360.
      if (flon.lt.-360.) flon = flon + 360.
      return
      end	
