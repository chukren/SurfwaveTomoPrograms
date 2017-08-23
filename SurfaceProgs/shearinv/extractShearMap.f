c  extract map of shear velocities at specified depth from 3Dcrust* models
c  after phase velocity inversion using shearz2.f
c
c  extractShearMap < extractinput
      character*70 shvelmodel, shveldpth(50)
      real*4 thick(50),beta(50), depth(50), dens(50)
      real*4 lon,lat
	read(*,*) shvelmodel
	read(*,*) ndepths
	do j = 1, ndepths
	  read(*,*) depth(j)	  
	  read(*,*) shveldpth(j)
	  open(10+j,file = shveldpth(j))
	enddo
c if density .gt. 3.2, interpolate velocities in mantle.  At density < 3.2, use 
c  actual velocities from layers so as to be able to represent crust
	
	open(10, file=shvelmodel)
	
	read(10,*) beglat,endlat,dlat,beglon,endlon,dlon
c	write(90,*) beglat,endlat,dlat,beglon,endlon,dlon
	nlat = (endlat-beglat)/dlat +1.01
        nlon = (endlon-beglon)/dlon + 1.01
        nxy = nlat*nlon

	do inxy = 1, nxy

	   read(10,*) lon, lat
           write(*,*) inxy, lon,lat
	   read(10,*) nlay
	   do i = 1, nlay
	     read(10,*) thick(i),dens(i),xx2,beta(i)
	   enddo
	   
	   do j = 1, ndepths
	     cumldepth = 0.0
	     do i = 1, nlay
	       cumldepth = cumldepth + thick(i)
	       if (depth(j).le.cumldepth) then
	         if ((dens(j).lt.3.2).or.(i.eq.1)) then
		   shvel = beta(i)
		 else
		   depthmid = cumldepth - thick(i)/2.0
		   if (depth(j).gt.depthmid) then
		     if (j.eq.ndepths) then
		       shvel = beta(i)
		     else
		       depthmid2 = cumldepth +thick(i+1)/2.0
		       shvel = (depth(j)-depthmid)*(beta(i+1)-beta(i))/
     1			 (depthmid2-depthmid) + beta(i)
     		     endif
		   else
		     if (dens(i-1).lt.3.2) then
		       shvel = beta(i) 
		     else
		       depthmid2 = cumldepth -thick(i)
     1                             -thick(i-1)/2.0
      		       shvel =(depth(j)-depthmid2)*(beta(i)-beta(i-1))/
     1                    (depthmid-depthmid2) + beta(i-1)
                     endif
		   endif
		 endif
		 write(10+j,*) lon,lat,shvel
		 go to 100
	       endif
	     enddo
100	     continue
	   enddo
	 enddo
	 do j = 1, ndepths
	   close(10+j)
	 enddo
	 close(10)
	 end
		   
	   	
