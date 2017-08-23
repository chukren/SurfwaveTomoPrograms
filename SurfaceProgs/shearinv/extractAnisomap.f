c Extracts cs2th and sin2th terms to plot vectors of anisotropy at specified depths
c using mapscriptSchanges
c
c extractAnisomap < extracAnisomapinp

       character*70 anisomodel, anisodepth(50)
       real*4 thick(50), cs2th(50), sin2th(50)
       real*4 lat, lon, depth(50),phi(50), amp1(50)
       real*4 slength(50)
       pi = 3.14159
       conv2deg = 180./pi
       scalefac = 4.0
       
       read(*,*) anisomodel
       read(*,*) ndepths
       do j = 1, ndepths
         read(*,*) depth(j)
	 read(*,*) anisodepth(j)
	 open(10+j, file = anisodepth(j))
       enddo
       
       open(10, file = anisomodel)
       
       read(10,*) beglat,endlat,dlat,beglon,endlon,dlon
       nlat = (endlat-beglat)/dlat +1.01
       nlon = (endlon-beglon)/dlon + 1.01
       nxy = nlat*nlon
       
       do inxy = 1, nxy
       
          read(10,*) lon, lat
	  read(10,*) nlay
	  
	  do i = 1, nlay
	     read(10,*) thick(i), cs2th(i), sin2th(i)
	     a1 = cs2th(i)
	     a2 = sin2th(i)
c phi is fast direction, theta is the azimuth from north	     
	     phi(i) = 0.5*atan2(a2,a1)
	     amp = sqrt(a1*a1+a2*a2)
	     amp1(i) = amp
	     
	     phi(i) = phi(i)*conv2deg
             slength(i) = scalefac*amp1(i)
	     write(*,*) conv2deg, phi(i), slength(i)
	  enddo
	  
	  do j = 1, ndepths
	     cumldepth = 0.0
	     do i = 1, nlay
	        cumldepth = cumldepth + thick(i)
		if (depth(j).le.cumldepth) then
		   depthmid = cumldepth - thick(i)/2.0
		   if (depth(j).gt.depthmid) then
		      if (j.eq.ndepths) then
		         phi2 = phi(i)
c			 amp2 = amp1(i)
                         slength2 = slength(i)
		      else
		         depthmid2 = cumldepth + thick(i+1)/2.0
			 phi2 = (depth(j)-depthmid)*(phi(i+1)-phi(i))/
     1                   (depthmid2-depthmid)+phi(i)
c                         amp2 = (depth(j)-depthmid)*(amp1(i+1)-amp1(i))/
c     1                   (depthmid2-depthmid)+amp1(i)
                          slength2 = (depth(j)-depthmid)*(slength(i+1)-
     1                    slength(i))/(depthmid2-depthmid)+slength(i)
                      endif
		   else
		   depthmid2 = cumldepth - thick(i)-thick(i-1)/2.0
		   phi2 = (depth(j)-depthmid2)*(phi(i)-phi(i-1))/
     1               (depthmid-depthmid2) + phi(i-1)
c                   amp2 = (depth(j)-depthmid2)*(amp1(i)-amp1(i-1))/
c     1               (depthmid-depthmid2) + amp1(i-1)
                   slength2 = (depth(j)-depthmid2)*(slength(i)-
     1               slength(i-1))/(depthmid-depthmid2) + slength(i-1)
		   endif
c		   slength2 = scalefac*amp2*10
		write(10+j,*) lon, lat, phi2, slength2
		go to 100
             endif
	     enddo
100          continue
            enddo
	   enddo
	   
	   do j = 1, ndepths
	      close(10+j)
	   enddo
	   close(10)
	   
	   end	     
           
		      
		
		   
		
	  
	  
