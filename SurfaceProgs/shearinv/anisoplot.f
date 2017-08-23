c Converts cs2th and sin2th terms to plot anisotropy on a vertical profile
c with perspective viewing
c
c anisoplot < anisoplotinp

       character*70 anisomodel, anisoout, anisoout2
       character*70 anisoout3
       real*4 thick(50), cs2th(50), sin2th(50)
       real*4 lat, lon, depth(50), phi(50), amp1(50)
       real*4 latin, lonin, slength(50),depth2(50)
       
       pi = 3.14159
       conv2deg = 180./pi
       
       read(*,*) anisomodel
       read(*,*) scalefac
       read(*,*) lonin, latin
       read(*,*) ndepths
       do j = 1, ndepths
          read(*,*) depth(j)
       enddo
       read(*,*) anisoout
       read(*,*) anisoout2
       read(*,*) anisoout3
       
       open(10, file = anisomodel)
       open(20, file = anisoout)
       open(30, file = anisoout2)
       open(40, file = anisoout3)
       
       write(20,*) '0 0 0 200 90 3.2'
       write(30,*) '0 0 0 200 90 3.2'
       write(40,*) '0 0 0 -1 90 3.2'
       
       read(10,*) beglat,endlat,dlat,beglon,endlon,dlon
       nlat = (endlat-beglat)/dlat +1.01
       nlon = (endlon-beglon)/dlon + 1.01
       nxy = nlat*nlon
       
       do inxy = 1, nxy
          read(10,*) lon, lat
	  if ((lon.eq.lonin).and.(lat.eq.latin)) then
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
	  enddo
	  
	 do j = 1, ndepths
	     cumldepth = 0.0
	     do i = 1, nlay
	        cumldepth = cumldepth + thick(i)
		if ((depth(j).le.cumldepth).and.(depth(j)
     1                     .gt.(cumldepth-thick(i)))) then
                   indlayer = i
		endif
	      enddo
	      depth2(j) = depth(j)/10
	      phi2 = phi(indlayer)
	      slength2 = slength(indlayer)
c phi3 is rotated 90 degrees for plotting purposes
              phi3 = phi2+90
	      amp3 = (slength2/scalefac)*100	      
	      write(20,*) '0 0' , depth2(j), phi2, phi2, slength2
	      write(30,*) '0 0' , depth2(j), phi2, phi3, slength2
	      write(40,*) '0 0' , depth2(j), amp3, phi2, slength2
	   enddo
	  
	  
	  else
	  read(10,*) nlay
	  do i = 1, nlay
	     read(10,*) thick(i), cs2th(i), sin2th(i)
	  enddo
	  endif
	  enddo
	  
	  close(10)
	  close(20)
	  
	  end
	     
	     
