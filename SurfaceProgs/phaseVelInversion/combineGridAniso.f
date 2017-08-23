c  prepare input for shear velocity inversion by combining phase velocities 
c  from inversions at individual periods into single file with velocities at 
c  all periods for each location.  ( added changes back to predicted
c  phase velocities from crustal model in files like gridphase.060.564.40.71)
c  
c  combineGridAniso < combineGridAnisoInp
      parameter (npermax = 40, ngridmax = 20000)
      character*70 foutput,finput(npermax)
      real*4 lon(ngridmax),lat(ngridmax)
      real*4 c2vel(npermax,ngridmax),s2vel(npermax,ngridmax)
      real*4 c2stddv(npermax,ngridmax),s2stddv(npermax,ngridmax)
      read(*,*) foutput
      read(*,*) nper
      do i = 1, nper
        read(*,*) finput(i+10)
	open(i+10, file = finput(i+10))
	read(i+10,*) beglat,endlat,dlat,beglon,endlon,dlon
      enddo
      nlat = (endlat-beglat)/dlat +1.01
      nlon = (endlon-beglon)/dlon + 1.01
      nxy = nlat*nlon
      do i = 1,nper
        do ii = 1, nxy
          read(i+10,*) lon(ii),lat(ii),
     1      c2vel(i,ii),c2stddv(i,ii),s2vel(i,ii),s2stddv(i,ii)
	enddo
      enddo  
      open (50,file = foutput) 
      write(50,*) nper   
      do ii = 1,nxy
        write(50,*) lon(ii),lat(ii)
        do i = 1,nper
	  write(50,*) c2vel(i,ii),c2stddv(i,ii),
     1                s2vel(i,ii),s2stddv(i,ii)
	enddo
      enddo
      
      do i = 1,nper
        close(i+10)
      enddo
      
      close(50)
      
      end
      
