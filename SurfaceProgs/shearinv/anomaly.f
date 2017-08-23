C Calculation the shear velocity anomaly
C Read in  shearvel*km.dat files and write out anomaly*.grd
	
	parameter (nlay = 12, nxy = 4686)
        character*8 depths(nlay)
        real*4 lon,lat,vel,ave,aver(nlay),anomaly,stdde,percentage
	

	depths(5) = '60.00000'
      	depths(6) = '80.00000'
	depths(7) = '100.0000'
	depths(8) = '120.0000'
	depths(9) = '140.0000'
	depths(10) = '160.0000'
      	depths(11) = '180.0000'
	depths(12) = '200.0000'
c	depths(12) = '220.0000'
c	depths(13) = '240.0000'
c	depths(14) = '265.0000'
c	depths(15) = '295.0000'
	

c	aver(5) = 4.3365
c	aver(6) = 4.3061
c	aver(7) = 4.1619
c       aver(8) = 4.0739
c	aver(9) = 4.0793
c	aver(10) = 4.1483
c	aver(11) = 4.2273
c	aver(12) = 4.3124
c	aver(12) = 4.4054
c	aver(13) = 4.4543
c	aver(14) = 4.5087
c	aver(15) = 4.5791
	
	
        Do i = 5, nlay

           open(1, file = 'shearvel.'//depths(i)//'km.dat')
	   open(2, file = 'shear_anomaly.'//depths(i)//'km.dat')
           open(3, file = ' stdde.'//depths(i)//'km.dat')
           open(4, file = ' stdde_percentage.'//depths(i)//'km.dat')
           sum=0.0
           Do k = 1, nxy
              read(1,*) xx, xx, vel
              sum=sum+vel
           enddo
           ave=sum/nxy
           write(*,*) depths(i),ave

           rewind(1)
           Do k = 1, nxy
             read(1,*) lon, lat, vel
             read(3,*) lon, lat, stdde
             anomaly = (vel - ave)/ave*100
             percentage= stdde/ave*100
             write(2,*) lon, lat, anomaly
             write(4,*) lon, lat, percentage
           Enddo

           close(1)
           close(2)
        Enddo  
 

	END




