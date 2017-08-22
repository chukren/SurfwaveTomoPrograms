c  Generate predicted phase velocities for 3D structure.

c differs from phsevel.f in having number of lower mantle layers read in,
c  instead of presumed, so that can change depth of later inversion easier.

c  Forward Problem  
c  Study region (-120~-106 & 22~35)
c  Predict phase velocity on 2D grids base on known 3D structure.
c  Modified from Don's inversion program 'shearz.f' in 'saito' directory. 

c  Modify statements between '****' lines to adjust to different study. 
c  May not indicate all the places where changes are needed

c  The 3D structure (density, vp, vs) used here is output from gencruststrt2.f

c  Phase velocities are calculated by s_saito.f.

c  ./phsevel2 < phsevel2.input

********************************Modify here**************************
c  nlay: number of layers in the lithosphere and upper mantle 
      parameter (npermax=30,maxlay = 100)  
*********************************************************************                                
     
      real*4	depth(maxlay)
      real*4	zero,thick(maxlay),dens(maxlay),density
      real*4    depmant(maxlay),densmant(maxlay),almant(maxlay),betmant(maxlay)	
      real*4    psait(npermax),phasait(npermax),bandsait(npermax)
      real*4    depsait(npermax),betsait(npermax)
      real*4    lat,lon,lontemp, lattemp 
      real*4    avgvel(npermax)
      
      double precision ccpred(npermax)
      double precision perd(npermax),xxd1,xxd3,xxp
      double precision dca(maxlay,npermax),dcb(maxlay,npermax)
    
      integer*4 iter,icnt,irow,icol,indx(npermax)
      integer*4 nparam,i,j,k,kk,xxd2,xxd4
      
      character*70 xxd0, flowmantle
      character*70 fcrust,fpredfase,fgrid(npermax)
	 
	read(*,*) fcrust
	read(*,*) flowmantle
	read(*,*) fpredfase
	
	read(*,*) nper2, depth1, error1,ijk1,ijk2,ijk3,topd,bottomd
c  These are variables controlling the search for phase velocity (roots to equations) for
c  given structure.  nper2= # periods, depth1 default depth extent, error1, criterion for 
c  ending search, ijk1 = 4 is Rayleigh spherical earth, 3 is Love spherical earth, ijk2 = 0, ijk3 = 2
c  topd and bottomd control depth range in which partial derivatives are output  


        open (10, file = fcrust)
	open( 20, file = fpredfase)
	
	read(*,*) nper
	do i = 1, nper
	  read(*,*) fgrid(i)
	  avgvel(i) = 0.0
	enddo

********************************Modify here**************************
	read(10,*) beglat,endlat,dlat,beglon,endlon,dlon
	write(20,*) beglat,endlat,dlat,beglon,endlon,dlon
	nx = (endlon-beglon)/dlon + 1.01
	ny = (endlat-beglat)/dlat +1.01
         nxy =  nx * ny
********************************************************************* 
         DO 100 inxy = 1, nxy
 
c+++++++++++++++++++++Step 1++++++++++++++++++++++++++++
     
c  Build input file 'scalifornia_refer.dat' for 's_saito.f'. 
      open(30, file = 'scalifornia_refer.dat',status='unknown')	    
c  Read in upper 210 km model
        read(10,*) lontemp, lattemp
	write(*,*) lontemp,lattemp
        read(10,*) nlay
********************************Modify here**************************      
c  54 is twice the number of layers in the lower mantle (below 210 km)
c  check file 'lowmantle.d'
c  nlayall: 2X layers, including low mantle parts 
c      nlayall= 2*nlay + 54
*********************************************************************  
c  Read in lower mantle model (below 210 km or whatever) 	  
c  unit=40: lowmantle.d
c  20 to 21 line P velocity digitals increse by 1  
        open( 40, file = flowmantle)
	read(40,*) nlowlay
        do i=1,nlowlay                                                               
	   read(40,*) depmant(i),densmant(i),almant(i),betmant(i)            
	enddo
      nlayall= 2*nlay + nlowlay

	 write(30,*)'    0    0    0'
         write(30,*) '   1'
         write(30,*)'df              col # = 47'
         write(30,*) nlayall, '    6371.0                    982.0                          
     r     0 0 1'           
         write(30,*)'(F6.2,3F8.4)'
     
        zero = 0.
c	write(*,*) nlay
	do ii=1,nlay
           read(10,*)    thicktemp,dentemp,ptemp,stemp 
c	   write(*,*)   thicktemp,dentemp,ptemp,stemp    
	   write(30,210) zero,dentemp,ptemp,stemp
	   write(30,210) thicktemp,dentemp,ptemp,stemp 
	enddo


        do i=1,nlowlay                                                               
	   write(30,210) depmant(i),densmant(i),almant(i),betmant(i)            
	enddo
	
c  nlowlay is 2x the layers in the lower mantle (0 thickness to 
c  indicate discontinuity, followed by actual thickness of layer)
c	do i=21,54
c           read(40,*) depmant(i),densmant(i),almant(i),betmant(i)
c           write(30,210) depmant(i),densmant(i),almant(i),betmant(i)
c	enddo

	close(40)
        
	write(30,*) '   4   2872.3  0' 
	write(30,*) '   1'  

********************************Modify here**********************************		
c  18 is number of periods and 210 is the maximum depth where the data can constrain 
c	write(30,*) '  18  700.0    1.00E-06  4 0 2   0.0        210.0'  
c  CAUTION  important that output below has format(spacing) identical to line above
   	write(30,310) nper2, depth1, error1,ijk1,ijk2,ijk3,topd,bottomd
310     format(I5,F7.1,E12.2,1x,3I2,F6.1,F12.1)
 
c	write(30,*) '  18  700.0    1.00E-06  4 0 2   0.0        210.0'                  
	write(30,*) ' (6F10.6,2I2)'
*****************************************************************************
        open(41, file='perstart2.d')
        do i = 1, nper
        read(41,*) psait(i),phasait(i),bandsait(i),betsait(i),depsait(i)
        write(30,230) psait(i),phasait(i),bandsait(i),betsait(i),
     r   depsait(i)
        enddo
        close(41)
        close(30)
c+++++++++++++++++++++Step 2+++++++++++++++++++++++++++++++++++
c calculate predicted phase velocities from s_saito.f
        call s_saito

c+++++++++++++++++++++Step 3+++++++++++++++++++++++++++++++++++
c read in predicted data and write into 'gridphase.dat'
c here c(w) obtained from s_saito which also calculates 
c partial derivatives,dcb,dca

	 open(16, file='DERIV.DATA')
 
c         write(*,*) lontemp, lattemp
         write(20,240) lontemp, lattemp 
	 do i=1,nper

           read(16,*) xxd0                           !!!!unit=16:  DERIV.DATA
c	   write(*,*) xxd0
           read(16,*) xxd1
c	   write(*,*) xxd1
           read(16,*) xxd2,xxd3,perd(i),ccpred(i)
           read(16,*) xxd4
	   
           write(20,250) ccpred(i)
	   avgvel(i) = avgvel(i) + ccpred(i)
c           write(*,*) nlay, perd(i)           
 	   do j=1,nlay
	    read(16,*) depth(j),density,dca(j,i),dcb(j,i)
	   enddo
c           write (*,*) depth(1),perd(i)
	enddo
        close (16)

100   ENDDO
	
      close(10)
      close(20)

c200   format (1x,f5.2,2x,f6.4,2x,f6.4,2x,f6.4)
210   format (1x,f5.2,2x,f6.4,2x,f7.4,1x,f6.4)
230   format (3x,f7.3,5x,f5.3,5x,f5.3,5x,f5.3,3x,f8.3)
240   format (2f8.2)
250   format (3x,f10.4)

c Reassemble predicted phase velocity file to create grid phase file for each period
	open (20,file = fpredfase)
	read(20,*) beglat,endlat,dlat,beglon,endlon,dlon
	do i = 1,nper
	  ifn = 20 + i
	  open(ifn, file = fgrid(i), status = 'unknown')
	  write(ifn,*) beglat,endlat,dlat,beglon,endlon,dlon
	  avgvel(i) = avgvel(i)/nxy
	  write(ifn,*) avgvel(i)
	enddo
        DO i =1, nxy
           read(20,*) lon,lat
           DO j = 1, nper
	     ifn = 20 +j
             read(20,*) phc
	     write(ifn,*) lon,lat,phc
           ENDDO
        ENDDO 
	do i = 1, nper
	  ifn = 20 + i
	  close(ifn) 
	enddo
	close(20)
         
      end

      
      
