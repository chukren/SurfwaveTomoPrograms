c Extracts the midpoint of each layer of a lat/lon point to use
c in a vertical profile
c
c extractVertProf < extractVertProfinp

       character*70 shvelmodel, profout
       real*4 thick(50),beta(50), depth(50), dens(50)
       real*4 lon,lat, halfthick(50)
       
       
       read(*,*) shvelmodel
       read(*,*) prolon, prolat
       read(*,*) profout
       
       open(10, file = shvelmodel)
       open(20, file = profout)
       
       read(10,*) beglat,endlat,dlat,beglon,endlon,dlon
       nlat = (endlat-beglat)/dlat +1.01
       nlon = (endlon-beglon)/dlon + 1.01
       nxy = nlat*nlon
	
       do inxy = 1, nxy
          cumldepth = 0.0
          read(10,*) lon, lat
	  if ((lon.eq.prolon).and.(lat.eq.prolat)) then
	     read(10,*) nlay
	     do i = 1, nlay
	        read(10,*) thick(i),dens(i),xx2,beta(i)
             enddo
	     do i = 1, nlay
	        halfthick(i) = thick(i)/2
	        depth(i) = cumldepth + halfthick(i)
	        write(20,*) depth(i), beta(i)
		if (i.lt.nlay) then
		if ((dens(i).lt.3.2).and.(dens(i+1).gt.3.2)) then
		  depth2 = depth(i)+halfthick(i)
		  write(20,*) depth2, beta(i)
		  write(20,*) depth2, beta(i+1)
		endif
		endif
	        cumldepth = cumldepth+thick(i)
	     enddo
	  else
	     read(10,*) nlay
	     do i = 1, nlay
	        read(10,*) thick(i),dens(i),xx2,beta(i)
	     enddo
	  endif
       enddo
       
       end
	  
