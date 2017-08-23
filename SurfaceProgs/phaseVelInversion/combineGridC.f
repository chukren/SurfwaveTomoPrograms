c  prepare input for shear velocity inversion by combining phase velocities 
c  from inversions at individual periods into single file with velocities at 
c  all periods for each location.  ( added changes back to predicted
c  phase velocities from crustal model in files like gridphase.060.564.40.71)
c  
c  combineGridC < combineGridCinp
      parameter (npermax = 40, ngridmax = 20000)
      character*70 foutput,finput(npermax)
      real*4 lon(ngridmax),lat(ngridmax),vel(npermax,ngridmax)
      real*4 stddv(npermax,ngridmax)
      read(*,*) foutput
      read(*,*) nper
      do i = 1, nper
        read(*,*) finput(i+10)
	open(i+10, file = finput(i+10))
	read(i+10,*) beglat,endlat,dlat,beglon,endlon,dlon
	read(i+10,*) preunifvel
      enddo
      nlat = (endlat-beglat)/dlat +1.01
      nlon = (endlon-beglon)/dlon + 1.01
      nxy = nlat*nlon
      do i = 1,nper
        do ii = 1, nxy
          read(i+10,*) lon(ii),lat(ii),vel(i,ii),stddv(i,ii)
	enddo
      enddo  
      open (50,file = foutput)    
      do ii = 1,nxy
        write(50,*) lon(ii),lat(ii)
        do i = 1,nper
	  write(50,*) vel(i,ii),stddv(i,ii)
	enddo
      enddo
      
      do i = 1,nper
        close(i+10)
      enddo
      
      close(50)
      
      end
      
