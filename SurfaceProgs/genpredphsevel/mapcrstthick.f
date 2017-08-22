c  quick routine to extract crustal thickness from 3Dcrust starting models
c  mapcrstthick < mapcrstinp
      character*70 fcrust,fthick
      real*4 lontemp,lattemp
	 
	read(*,*) fcrust
       read(*,*) fthick
        open (10, file = fcrust)
	open( 20, file = fthick)

	read(10,*) beglat,endlat,dlat,beglon,endlon,dlon
	nx = (endlon-beglon)/dlon + 1.01
	ny = (endlat-beglat)/dlat +1.01
         nxy =  nx * ny
********************************************************************* 
      do inxy = 1, nxy
        read(10,*) lontemp, lattemp
        read(10,*) nlay
	icrust = 1
	cthick = 0.0
	do ii=1,nlay
           read(10,*)    thicktemp,dentemp,ptemp,stemp 
	   if ((icrust.eq.1).and.(stemp.lt.4.0)) then
	     cthick = cthick +thicktemp
	   else
	     icrust = 0
	   endif
	 enddo
         write (20,*) lontemp,lattemp,cthick
       enddo
       close (20)
       end
       
