c  generate crustal starting model for SWCal+B&R to be used for generating 
c  predicted phase velocities that will be the starting point for inversion for 
c  phase velocities from multi-plane wave method
c
c  Assume will generate on regular grid of lat/long points  (could read in
c   points to be used instead.
c
c  Same as gencruststrt2.f except replace Levander's moho model with ours and
c  use linear interpolation instead of nearest value

c  assume start in lower left (SW) corner

      real*4 beglat,endlat, latinc,beglon,endlon, loninc
      real*4 coastlat(90),coastlon(90)
      real*4 costlon(90,20),wdepth(90,20)
      real*4 thickl(50),density(50),vp(50),vs(50)
      real*4 defthick(50),defdens(50),defvp(50),defvs(50)
      real*4 shlinlat(50),shlinlon(50)
      real*4 sllat(29,45),sllon(29,45),slvp(29,45,8),slvpvs(29,45,8)
      real*4 moholat(21000),moholon(21000),mdepth(21000)
      real*4 moslat(4000),moslon(4000),mosrayl8(4000)
      integer*4 nnc(50)
      character*70 coastfile, outputfile, mantlefile,linshearfile
      character*70  mohofile, defaultfile, rayl8secfile
      common /moschetti1/ moslat,moslon,mosrayl8,nmos
            
      read(*,*) beglat,endlat,latinc, beglon,endlon, loninc
      read(*,*) coastfile
      read(*,*) outputfile
      read(*,*) mantlefile
      read(*,*) linshearfile
      read(*,*) mohofile
      read(*,*) defaultfile
      read(*,*) rayl8secfile
      
      open(10, file = coastfile)
      open(11, file = outputfile)
      open(12, file = mantlefile)
      open(13, file = linshearfile)
      open(14, file = mohofile)
      open(15, file = defaultfile)
      open(16, file = rayl8secfile)
      open (17, file = 'testmoho')
      
      write(11,*) beglat,endlat,latinc, beglon,endlon, loninc
      
c read in Moshcetti etal 8 s period Rayleigh wave solution 
      read(16,*) nmos
      do i = 1,nmos
        read(16,*) moslon(i),moslat(i),mosrayl8(i)
	moslon(i) = moslon(i) - 360.
      enddo     

c read in bounds of Lin,Shearer et al.(2007) model that is reasonably constrained
      read(13,*)  nshlin
      do i = 1, nshlin
        read(13,*) shlinlon(i),shlinlat(i)
      enddo
c  read in Lin, Shearer model (15 km grid)
      do k = 1,8
       do j = 1, 45
        do i = 1,29
	  read(13,*)  sllon(i,j),sllat(i,j),dummy,dummy2,dummy3,
     1        slvp(i,j,k),slvpvs(i,j,k)
        enddo
       enddo
      enddo
c read in default crust/mantle model
      read(15,*) ndef
      do i = 1, ndef
        read(15,*) defthick(i),defdens(i),defvp(i),defvs(i)
      enddo      
c  read in Levander's moho depth model
c      nmoho = 20053
c      do i = 1, nmoho
c        read(14,*) moholat(i),moholon(i),mdepth(i)
c      enddo
c  read in our crustal model from receiver functions (highly edited 
c  Crotwell and Owens IRIS library)
      nmoho = 0
      do i =1, 100000
        read(14,*) moholon(i),moholat(i),mdepth(i)
	if (moholat(i).gt.91.) go to 200
	nmoho = nmoho + 1
      enddo
200   continue

c  read default mantle structure 
      read(12,*) nmantle
      do i = 1, nmantle
        read(12,*) thickl(i),density(i),vp(i),vs(i)
      enddo
      
      read(10,*) ncoast
c ncoast is number of coastal points and number of EW profiles on which water
c  depth is described.  
      do i = 1, ncoast 
        read(10,*) coastlat(i),coastlon(i)
      enddo
      do i = 1, ncoast
        read (10,*) nnc(i)
	do j = 1, nnc(i)
          read (10,*) dummy, costlon(i,j), wdepth(i,j)
c	  write(*,*) dummy, costlon(i,j), wdepth(i,j)
	enddo
      enddo  
      do slat = beglat,endlat, latinc
        do slon = beglon,endlon, loninc
	  nregion = 0
c  need to find which model regime point is in.  
c
c  first, look for offshore (100 m water or deeper), from crude outline
	  if ((slat.gt.coastlat(ncoast)).or.(slat.lt.coastlat(1))) then
	    write(*,*) slat,slon, 'slat out of coastal range'
	  endif
c  test for location between profiles and then for location west of coast
	  do ic = 1, ncoast-1
	    if ((slat.ge.coastlat(ic)).and.(slat.lt.coastlat(ic+1))) then 
	      coint = (slat-coastlat(ic))/(coastlat(ic+1)-coastlat(ic))
	      colonint = coint*coastlon(ic+1)+ (1.0-coint)*coastlon(ic)
c	      write (*,*) slat,slon,ic, coint, colonint
	      if (slon.le.colonint) then
	        nregion = 1
	        if (slon.lt.costlon(ic,1)) wd1= wdepth(ic,1)
		if (slon.lt.costlon(ic+1,1)) wd2 = wdepth(ic+1,1)
		if (slon.gt.costlon(ic,nnc(ic))) wd1= wdepth(ic,nnc(ic))
		if (slon.gt.costlon(ic+1,nnc(ic+1)))
     1		    wd2=wdepth(ic+1,nnc(ic+1))
     		do jj = 1, nnc(ic)-1
		  if ((slon.ge.costlon(ic,jj)).and.
     1		      (slon.le.costlon(ic,jj+1))) then
		    wd1int = (slon-costlon(ic,jj))/
     1		      (costlon(ic,jj+1)-costlon(ic,jj))
		    wd1 = wd1int*wdepth(ic,jj+1) 
     1                  +(1.0-wd1int)*wdepth(ic,jj)
     		  endif
		  if ((slon.ge.costlon(ic+1,jj)).and.
     1		      (slon.le.costlon(ic+1,jj+1))) then
		    wd2int = (slon-costlon(ic+1,jj))/
     1		      (costlon(ic+1,jj+1)-costlon(ic+1,jj))
		    wd2 = wd2int*wdepth(ic+1,jj+1) 
     1                  +(1.0-wd2int)*wdepth(ic+1,jj)
     		  endif
		enddo
		waterd = coint*wd2 +(1.0-coint)*wd1
		call wdcrustgen(waterd,slat,slon,thickl,density,vp,vs,
     1              zmoho,nmantle)
	      endif
	    endif
	  enddo
	  if (nregion.eq.0) then
c  If location not off the coast, estimate moho thickness from receiver function
c  values
c  Instead of nearest point for Levander, in this version,
c  use weighted average of points within 0.75 or if none within 0.75 deg, within 1.0 deg
c  weights range from 1.0 at site to 0 at limit distance, inversely proportional to 
c  distance
	    wgtsum = 0.0
	    thicksum = 0.0
	    do i = 1,nmoho
	      dist2 = (slat-moholat(i))**2 + (slon-moholon(i))**2
	      dist = sqrt(dist2)
	      if (dist.lt.0.75) then
	        wgti = 0.75-dist
	        wgtsum = wgtsum + wgti
		thicksum = thicksum + mdepth(i)*wgti
              endif
	    enddo
	    if (wgtsum.ne.0.0) then
	      zmoho = thicksum/wgtsum
	    else
	      do i = 1,nmoho
	        dist2 = (slat-moholat(i))**2 + (slon-moholon(i))**2
	        dist = sqrt(dist2)
	        if (dist.lt.1.0) then
		  wgti = 1.0-dist
	          wgtsum = wgtsum + wgti
		  thicksum = thicksum + mdepth(i)*wgti
                endif
              enddo
	      if (wgtsum.ne.0.0) then
	        zmoho = thicksum/wgtsum
	      else
c  If too far from any receiver-function-constrained point, use 31 km crust.
	        zmoho = 31.0
	      endif
	    endif
	        
c  Test if within Lin, Shearer etal. model region
c  Uses mercator x,y test, not spherical earth region with great circle bounds
  	    inlin = 0
	    call polygin(shlinlon,shlinlat,nshlin,slon,slat,inlin)
	    if (inlin.eq.1) then
	      nregion = 2
	      write (*,*) slat,slon,' shearer'
c  find nearest point in shearer, adopt values - don't worry about interpolation since
c  15-km grid is finer than our resolution
              distmin = 1.0E+30
  	      do i = 1,29
	        do j = 1,45
	          dist2 = (slat-sllat(i,j))**2 + (slon-sllon(i,j))**2
	          dist = sqrt(dist2)
                  if (dist.lt.distmin) then
		    distmin = dist
		    imin = i
		    jmin = j
		  endif
	        enddo
	      enddo
	      call socalgen(slat,slon,zmoho,imin,jmin,slvp,slvpvs,
     1         thickl,density,vp,vs,nmantle)
	    endif
c  If point not within any other region, use default crustal model (equivalent
c  to Lin,Shearer except 6.8 km/s lower crust) and receiver function crustal thickness
	    if (nregion.eq.0) then
	      call defaultgen(slat,slon,zmoho,defthick,
     1	        defdens,defvp,defvs,ndef)
	    endif
	  endif
	  write(17,300) slon,slat,zmoho
300   format (2F9.3,1x,f4.1)
	enddo
      enddo
      close(unit=10)
      close(unit=11)
      close(unit=12)
      close(unit=17)
      stop
      end
      
      subroutine defaultgen(slat,slon,zmoho,defthick,
     1	        dfdens,dfvp,dfvs,ndef)
      real*4 defthick(50),defdens(50),defvp(50),defvs(50),dpth(50)
      real*4 dfdens(50),dfvs(50),dfvp(50)
      real*4 drho(3),dvp(3),dvs(3)
c  find layer that moho lies in - make 1 km the minimum layer thickness
c  This routine specific to  default crust with 6 layers
      do i = 1,ndef
        defdens(i) = dfdens(i)
	defvs(i) = dfvs(i)
	defvp(i) = dfvp(i)
      enddo
      depth = 0.0
      write(11,101) slon,slat
      write(11,*) ndef
c  modify crustal properties upper 3 layers based on 8-sec Rayleigh wave
c  dispersion data from Moschetti et al
      call moschetti(slat,slon,drho,dvp,dvs)
      do i = 1,3
        defdens(i) = defdens(i) + drho(i)
	defvs(i) = defvs(i) + dvs(i)
	defvp(i) = defvp(i) + dvp(i)
      enddo      
      do i = 1, ndef
        depth = depth + defthick(i)
        dpth(i) = depth
      enddo
      do i = 2,ndef
	if ((zmoho.gt.dpth(i-1)).and.(zmoho.le.dpth(i))) then
           idepth = i
	endif
      enddo
c  if crust lies within original 6 layers, truncate crust and extend first
c  mantle layer properties upward , i.e. layer 7)
      if (idepth.le.6) then
        do i = 1, idepth-1
          write(11,100) defthick(i),defdens(i),defvp(i),defvs(i)
        enddo
        thick1 = zmoho - dpth(idepth-1)
        thick2 = dpth(idepth) - zmoho + defthick(idepth+1)
	if (thick1.lt.1.0) then
	  thick1 = 1.0
	  thick2 = defthick(idepth) -1.0 + defthick(idepth+1)
	endif
	write(11,100) thick1,defdens(idepth),defvp(idepth),
     1       defvs(idepth)
	write(11,100) thick2,defdens(7),defvp(7),defvs(7)
	if (idepth.le.4) then
	  do i = idepth+2,6
	    write(11,100) defthick(i),defdens(7),defvp(7),defvs(7)
	  enddo
	endif
	iminus = 0
	if (idepth.eq.6) iminus = 1
      endif
      if (idepth.gt.6) then
        do i = 1, 6
          write(11,100) defthick(i),defdens(i),defvp(i),defvs(i)
        enddo
	if (idepth.gt.7) then
          do i = 7, idepth-1
            write(11,100) defthick(i),defdens(6),defvp(6),defvs(6)
          enddo
	  thick1= zmoho - dpth(idepth-1)
	  thick2 = dpth(idepth) - zmoho + defthick(idepth+1)
	  if (thick1.lt.2.0) then
	    thick1 = 2.0
	    thick2 = defthick(idepth) -2.0 + defthick(idepth+1)
	  endif
	  write(11,100) thick1,defdens(6),defvp(6),
     1       defvs(6)
	  write(11,100) thick2,defdens(idepth+1),defvp(idepth+1),
     1       defvs(idepth+1)
          iminus = idepth - 5
	endif
	if (idepth.eq.7) then
	  thick1= zmoho - dpth(idepth-1)
	  thick2 = dpth(idepth) - zmoho + defthick(idepth+1)
	  if (thick1.lt.2.0) then
	    thick1 = 2.0
	    thick2 = defthick(idepth) -2.0 + defthick(idepth+1)
	  endif
	  write(11,100) thick1,defdens(6),defvp(6),
     1       defvs(6)
	  write(11,100) thick2,defdens(idepth+1),defvp(idepth+1),
     1       defvs(idepth+1)
	  iminus = 2
	endif	  
      endif      
      do i = 7+iminus, ndef
	write(11,100) defthick(i),defdens(i),defvp(i),defvs(i)
      enddo
100   format(F7.3,2x, 3(F6.3,3x))
101   format(F9.3,2x,F7.3)
      return
      end        

      subroutine moschetti(slat,slon,drho,dvp,dvs)
      real*4 moslat(4000),moslon(4000),mosrayl8(4000)
      real*4 drho(3),dvp(3),dvs(3)
      common /moschetti1/ moslat,moslon,mosrayl8,nmos
      
c  find points within 0.6 deg of slat,slon and average them.  Moschetti on 0.5 deg
c  grid.  Averaging provides both interpolation and smoothing - seems to be some
c  oscillation in their solution with largest values adjacent to smallest, so do a
c  little smoothing.
c 
c  background value for phase velocity in their model is 3.0772 km/s at 8 s period 
      weight = 0.0
      sum = 0.0
      do i = 1,nmos
        dist = sqrt((slat-moslat(i))**2 + (slon-moslon(i))**2)
	if (dist.le.0.6) then
	  wgt = 1.0/(dist+0.3) 
	  weight = weight + wgt
	  sum = sum + mosrayl8(i)*wgt
	endif
      enddo
      davg = sum/weight - 3.0772
c  limit variability of model
      if (davg.gt.0.2) davg = 0.2
      if (davg.lt.(-0.4)) davg = -0.4
c  get about 0.1 km/s change in phase velocity for 0.4 km/s change in Vs in upper
c  3 km, 0.2 km/s in next 3 km and 0.1 in next 4.  Not really linear, but assume
c  it is, with the limitations of linearity assumption ameliorated by the
c  limitation on variability in the model.  (This assumes dVp = dVs and
c  drho = dVp/2.2
       dvs(1) = 4.0*davg
       dvs(2) = 2.0*davg
       dvs(3) = 1.0*davg
       do i = 1,3
         dvp(i) = dvs(i)
	 drho(i)= dvs(i)/2.2
       enddo
       return
       end

      subroutine socalgen(slat,slon,zmoho,imin,jmin,slvp,slvpvs,
     1  thickl,density,vp,vs,nmantle)
      real*4 slvp(29,45,8),slvpvs(29,45,8)
      real*4 thick(5),cvp(5),cvs(5),cden(5)
      real*4 thickl(50),density(50),vp(50),vs(50)

c generate crustal layers and mantle layers down to 210 km - crustal layering
c  based on Lin,Shearer etal. tomo and Levanders's moho - usedfixed lower crust 
c  velocity below 22 km, since is not constrained anyway by tomo
c
c  crustal layers are 0-3, 3-6, 6-10, 10-15, 15-22, 22-moho
c  assumes mantle layers 10 km thick down to 90 km
      denlow = 3.05
      pvlow = 6.8
      svlow = 3.89
      thick(1) = 3.0
      thick(2) = 3.0
      thick(3) = 4.0
      thick(4) = 5.0
      if (zmoho.ge.22.0) then 
        nlayer = nmantle + 4
	nst = 5
	write(11,101) slon,slat
	write(11,*) nlayer
	thick(5) = 7.0
	do i = 1,5
	  cvp(i) = (slvp(imin,jmin,i)+slvp(imin,jmin,i+1))/2.0
	  avvpvs = (slvpvs(imin,jmin,i)+slvpvs(imin,jmin,i+1))/2.0
	  if (avvpvs.lt.1.65) avvpvs = 1.65
	  if (avvpvs.gt.2.00) avvpvs = 2.00
	  cvs(i) = cvp(i)/avvpvs
	  cden(i) = cvp(i)/2.2
	  if (cden(i).lt.2.0) cden(i) = 2.0
	  write(11,100) thick(i),cden(i),cvp(i),cvs(i)
	enddo
	if ((zmoho.gt.35.0).and.(zmoho.lt.45.0)) then
	  thick1 = (zmoho-22.0)/2.0
	  write(11,100) thick1,denlow,pvlow,svlow
	  write(11,100) thick1,denlow,pvlow,svlow
	  thick2 = 50.0-zmoho
	  write(11,100) thick2,density(4),vp(4),vs(4)
	endif
	if (zmoho.le.35.0) then
	  thick1 = zmoho - 22.0
	  thick2 = 40.0-zmoho
	  nst = 4
	  if (thick1.lt.1.0) then
	    thick1 = 1.0
	    thick2 = 40.0-23.0
	  endif
	  write(11,100) thick1,denlow,pvlow,svlow
	  write(11,100) thick2,density(4),vp(4),vs(4)
        endif
	if (zmoho.ge.45.0) then
	  thick1 = (zmoho-22.0)/2.0
	  write(11,100) thick1,denlow,pvlow,svlow
	  write(11,100) thick1,denlow,pvlow,svlow
	  thick2 = (60.0-zmoho)/2.0
	  write(11,100) thick2,density(4),vp(4),vs(4)
	  write(11,100) thick2,density(4),vp(4),vs(4)
	  nst = 6	  
        endif
	if (zmoho.gt.50.0) write(*,*) slat,slon, zmoho,
     1       ' crust over 50 km'	  	  
      endif
      if ((zmoho.ge.16.0).and.(zmoho.lt.22.0)) then
        nlayer = nmantle + 4
	nst = 3
	write(11,101) slon,slat
	write(11,*) nlayer
	thick(5) = zmoho-15.0
	thick1 = 30.0 - zmoho
	do i = 1,5
	  cvp(i) = (slvp(imin,jmin,i)+slvp(imin,jmin,i+1))/2.0
	  avvpvs = (slvpvs(imin,jmin,i)+slvpvs(imin,jmin,i+1))/2.0
	  if (avvpvs.lt.1.65) avvpvs = 1.65
	  if (avvpvs.gt.2.00) avvpvs = 2.00
	  cvs(i) = cvp(i)/avvpvs
	  cden(i) = cvp(i)/2.2
	  if (cden(i).lt.2.0) cden(i) = 2.0
	  write(11,100) thick(i),cden(i),cvp(i),cvs(i)
	enddo
	write(11,100) thick1,density(2),vp(2),vs(2)
      endif
      if (zmoho.lt.16.0) write(*,*) slat,slon, zmoho,
     1       ' crust under 16 km'
      do i = nst, nmantle
	write(11,100) thickl(i),density(i),vp(i),vs(i)
      enddo
100   format(F7.3,2x, 3(F6.3,3x))
101   format(F9.3,2x,F7.3)
      return
      end
	  	  
     
	
      
     
      subroutine wdcrustgen(waterd,slat,slon,thickl,density,vp,vs,
     1   zmoho,nmantle)

      real*4 thickl(50),density(50),vp(50),vs(50)
      
c generate crustal layers and mantle layers down to 210 km - crustal layering
c  based on waterdepth alone
c  Water layer directly interpolated value
c  If water depth = > than 3000 m, then use standard oceanic crust with 2 km
c  upper crust P velocity 5 km/s and 5 km lower crust velocity 6.8 km/s
c  Interpolate shallower depths to at 100m water,  constant uppermost crust,
c  6 km/s layer varying from 0 to 13 km thickness, and 6.8 km/s layer varying
c  from 5 to 13 km thickness.
c
c  assumes mantle layers 10 km thick down to 90 km (beginning at 10 km)
c  if crust thickens, thin or eliminate mantle layers to compensate 
      waterden = 1.03
      watersv = 0.00
      waterpv = 1.50
      thickup = 2.0
      denup = 2.45
      denmid = 2.80
      denlow = 3.05
      svup = 2.63
      svmid = 3.43
      svlow = 3.89
      pvup = 5.0
      pvmid = 6.0
      pvlow = 6.8
      waterdp = waterd/1000.
100   format(F7.3,2x, 3(F6.3,3x))
101   format(F9.3,2x,F7.3)
      if (waterd.ge.3000.) then
      	 thickm1 = thickl(1) - (waterd-3000.)/1000.
	 thicklow = 5.0
	 nlayer = nmantle + 3
	 nst = 2
	 zmoho = 7. + waterdp
	 write(11,101) slon, slat
	 write(11,*) nlayer
	 write(11,100) waterdp, waterden,waterpv,watersv
	 write(11,100) thickup, denup,pvup,svup
	 write(11,100) thicklow, denlow,pvlow,svlow
	 write(11,100) thickm1,density(1),vp(1),vs(1)
	 do i = nst, nmantle
	   write(11,100) thickl(i),density(i),vp(i),vs(i)
	 enddo
      endif
      if (waterd.lt.3000.) then
        waterint = (3000.-waterd)/2900.
        totcrust = 2.0 + waterint*(26.0-5.0) +5.0
	zmoho = totcrust+waterdp
	thickmid = waterint*13.0
	thicklow = waterint*(13.0-5.0) +5.0
	if (thickmid.lt.2.0) then
	  thickmid = 0.0
	  thicklow = totcrust - 2.0
	  thickm1 = 20.0 - waterdp -totcrust
	  nst = 2
	  nlayer = nmantle + 3
	  ifirst = 1	    
	  write(11,101) slon, slat
	  write(11,*) nlayer
	  write(11,100) waterdp, waterden,waterpv,watersv
	  write(11,100) thickup, denup,pvup,svup
	  write(11,100) thicklow, denlow,pvlow,svlow
	  write(11,100) thickm1,density(ifirst),vp(ifirst),vs(ifirst)
	  do i = nst, nmantle
	    write(11,100) thickl(i),density(i),vp(i),vs(i)
	  enddo
	endif
	if ((thickmid.ge.2.0).and.(thickmid.lt.7.0)) then
	  thickm1 = 20.0 - waterdp -totcrust
	  nst = 2
	  nlayer = nmantle + 4
	  ifirst = 1
	  if (thickm1.lt.3.0) then
	    thickm1 = thickm1 + 10.0
	    nlayer = nlayer -1
	    ifirst = 2
	    nst = 3
	  endif	    
	  write(11,101) slon, slat
	  write(11,*) nlayer
	  write(11,100) waterdp, waterden,waterpv,watersv
	  write(11,100) thickup, denup,pvup,svup
	  write(11,100) thickmid,denmid,pvmid,svmid
	  write(11,100) thicklow, denlow,pvlow,svlow
	  write(11,100) thickm1,density(ifirst),vp(ifirst),vs(ifirst)
	  do i = nst, nmantle
	    write(11,100) thickl(i),density(i),vp(i),vs(i)
	  enddo
	endif
	if (thickmid.ge.7.0) then
	  thickm1 = 30.0 - waterdp -totcrust
	  nst = 3
	  nlayer = nmantle + 5
	  ifirst = 2
	  thickmid2 = thickmid/2.0
	  thicklow2 = thicklow/2.0
	  if (thickm1.lt.4.0) then
	    thickm1 = thickm1 + 10.0
	    nlayer = nlayer -1
	    ifirst = 3
	    nst = 4
	  endif	    
	  write(11,101) slon, slat
	  write(11,*) nlayer
	  write(11,100) waterdp, waterden,waterpv,watersv
	  write(11,100) thickup, denup,pvup,svup
	  write(11,100) thickmid2,denmid,pvmid,svmid
	  write(11,100) thickmid2,denmid,pvmid,svmid
	  write(11,100) thicklow2, denlow,pvlow,svlow
	  write(11,100) thicklow2, denlow,pvlow,svlow
	  write(11,100) thickm1,density(ifirst),vp(ifirst),vs(ifirst)
	  do i = nst, nmantle
	    write(11,100) thickl(i),density(i),vp(i),vs(i)
	  enddo
	endif
      endif
      return
      end	 
      subroutine polygin(x,y,n,xtest,ytest,in)
c test whether point (xtest,ytest) lies within a polygon
c  n (x,y) points must be supplied sequentially outlining polygon
c
c  in = 0 outside, in = 1, inside
c
c  find max and min y values of each line segment and of total set
c  this part could be spun off and supplied to routine so it would
c  only have to be done once, speeding up process
      real*4 x(100),y(100),xmax(100),ymax(100),xmin(100),ymin(100)
      in =0
      xmaxx = -1.0E+30
      ymaxx = -1.0E+30
      xminn = 1.0E+30
      yminn = 1.0E+30
      do i = 1, n-1
        xmax(i) = amax1(x(i),x(i+1))
	ymax(i) = amax1(y(i),y(i+1))
	xmin(i) = amin1(x(i),x(i+1))
	ymin(i) = amin1(y(i),y(i+1))
	xmaxx = amax1(xmaxx,xmax(i))
	ymaxx = amax1(ymaxx,ymax(i))
	xminn = amin1(xminn,xmin(i))
	yminn = amin1(yminn,ymin(i))
      enddo
        xmax(n) = amax1(x(1),x(n))
	ymax(n) = amax1(y(1),y(n))
	xmin(n) = amin1(x(1),x(n))
	ymin(n) = amin1(y(1),y(n))
      if ((xtest.ge.xminn).and.(xtest.le.xmaxx)) then
        if ((ytest.ge.yminn).and.(ytest.le.ymaxx)) then
c  shoot ray in positive y direction from (xtest,ytest).  If intersects
c  even number (including zero) of bounding segments, point is outside
          inum = 0
	  do i = 1, n-1
	    if (x(i+1).ne.x(i)) then
	      yint = y(i)+(y(i+1)-y(i))*(xtest-x(i))/(x(i+1)-x(i))
	      if ((yint.ge.ymin(i)).and.(yint.le.ymax(i)).and.
     1             (yint.ge.ytest)) then
                inum=inum+1
	      endif
	    else
	      if ((xtest.eq.x(i)).and.((ytest.ge.ymin(i)).and.
     1               (ytest.le.ymax(i)))) then
     		inum = inum+1
	      endif
	    endif
	  enddo
	    if (x(1).ne.x(n)) then
	      yint = y(n)+(y(1)-y(n))*(xtest-x(n))/(x(1)-x(n))
	      if ((yint.ge.ymin(n)).and.(yint.le.ymax(n)).and.
     1             (yint.ge.ytest)) then
                inum=inum+1
	      endif
	    else
	      if ((xtest.eq.x(n)).and.((ytest.ge.ymin(n)).and.
     1               (ytest.le.ymax(n)))) then
     		inum = inum+1
	      endif
	    endif	  
	  diff = float(inum)/2.0 - int(inum/2)
	  if (diff.gt.0.00001) in = 1
	endif
      endif
      return
      end
       
