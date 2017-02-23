!==============================================================================  
!  generate crustal starting model for SWCal+B&R to be used for generating 
!  predicted phase velocities that will be the starting point for inversion for 
!  phase velocities from multi-plane wave method
!
!  Assume will generate on regular grid of lat/long points  (could read in
!  points to be used instead.
!
!  Same as gencruststrt2.f except replace Levander's moho model with ours and
!  use linear interpolation instead of nearest value

!  assume start in lower left (SW) corner
!
!  10/29/14  chukren  Princeton 
!       Modified to F90 format for Cascadia Initiative Project
!==============================================================================
module moschettil
  real :: moslat(4000),moslon(4000),mosrayl8(4000)
  integer :: nmos
end

module region
  integer :: nregion
  integer :: couldbe_seamnt
end
     
program gencrust
  use globalage
  use moschettil
  use region
  use topopara
  use slabpara
  use wUSpara
  use thermodel
  integer, parameter :: nlaymax = 50
  integer, parameter :: npcoast = 200
  integer, parameter :: nmaxm   = 21000
  integer, parameter :: nlinmax = 50 
  real, parameter :: slabmohomax = 38.0 

  real :: beglat,endlat,latinc,beglon,endlon,loninc
  real :: coastlat(npcoast),coastlon(npcoast)
  !real :: costlon(npcoast,20),wdepth(npcoast,20)
  real :: thickl(nlaymax),density(nlaymax),vp(nlaymax),vs(nlaymax)
  real :: defthick(nlaymax),defdens(nlaymax),defvp(nlaymax),defvs(nlaymax)

  !=== Lin et al. model ===
  real :: shlinlat(nlinmax),shlinlon(nlinmax)
  real :: sllat(29,45),sllon(29,45),slvp(29,45,8),slvpvs(29,45,8)
 
  !=== slab ===
  real :: slabdepth, ddslab, watershelf, zslab
 ! real :: waterd, sedthck, sdep

  !=== sediments ===
  real :: sed_vs

  real :: moholat(nmaxm),moholon(nmaxm),mdepth(nmaxm)
  integer :: nnc(50)
  integer :: ndef
  integer :: n_lat, n_lon,ii,jj

  character(len=70) :: coastfile, outputfile, mantlefile,SnHmantlefile
  character(len=70) :: mohofile, defaultfile, rayl8secfile
        
  read(*,*) beglat,endlat,latinc,beglon,endlon,loninc
  read(*,*) coastfile
  read(*,*) outputfile
  read(*,*) mantlefile
  read(*,*) SnHmantlefile
  read(*,*) mohofile
  read(*,*) defaultfile
  read(*,*) rayl8secfile
  
  open(20, file = coastfile)
  open(12, file = mantlefile)
  open(13, file = SnHmantlefile) ! not use
  open(14, file = mohofile)
  open(15, file = defaultfile)
  open(16, file = rayl8secfile)

  open(11, file = 'tempoutput')
  open(17, file = 'testmoho.dat')
  open(18, file = 'testregion.dat')
  open(19, file = 'testslabtopo.dat')
  open(23, file = 'testbasement.dat')
  open(25, file = 'testwater.dat')
  open(24, file = 'testsedvs.dat')
  
  !=== 
  !  read in Moshcetti etal 8 s period Rayleigh wave solution 
  !===
  write(*,*) 'start reading',rayl8secfile
  read(16,*) nmos
  do i = 1,nmos
    read(16,*) moslon(i),moslat(i),mosrayl8(i)
    moslon(i) = moslon(i) - 360.
  enddo     
  close(16)
  write(*,*) 'moshcetti done'

  !=== 
  !  read in default crust+mantle continental model
  !     defaultcrustCascade.dat
  !===
  write(*,*) 'start reading',defaultfile
  read(15,*) ndef
  do i = 1, ndef
    read(15,*) defthick(i),defdens(i),defvp(i),defvs(i)
  enddo  
  close(15)
  write(*,*) 'default continental model done!'  
  
  !===========================================================
  !  read in Levander's moho depth model
  !      nmoho = 20053
  !      do i = 1, nmoho
  !        read(14,*) moholat(i),moholon(i),mdepth(i)
  !      enddo
  !===========================================================

  !  read in our crustal model from receiver functions (highly edited 
  !  Crotwell and Owens IRIS library)
  write(*,*) 'start reading',mohofile
  nmoho = 0
  do i =1, 10000
    read(14,*) moholon(i),moholat(i),mdepth(i)
    if (moholat(i).gt.91.) go to 200
    nmoho = nmoho + 1
  enddo
200   continue
  write(*,*) 'continental moho done!'  
  

  !===  
  !  read default oceanic mantle structure 
  !  mantlemodel.dat
  !===
  write(*,*) 'start reading',mantlefile
  read(12,*) nmantle
  do i = 1, nmantle
    read(12,*) thickl(i),density(i),vp(i),vs(i)
  enddo
  close(12)
  write(*,*) 'oceanic mantle done!'  
  
  !===
  ! ncoast is number of coastal points 
  !=== 
  read(20,*) ncoast
  do i = 1, ncoast 
    read(20,*) coastlat(i),coastlon(i)
  enddo
  close(20)
  write(*,*) 'coast done!'  

  !=== 
  ! Read the water depth table: JdFtopo.bin  
  !===
  call readtopobin()
  write(*,*)'read water depth done!'

  !===
  ! Read the NOAA's sediment thickness table: JdFsed.xyz
  !===
  call read_sed_xyz()
  write(*,*) 'read sediment thickness done!'

  !===
  ! Read the global seafloor age database
  !===
  call readglobalcrustage()
  write(*,*) 'read seafloor age done!'

  !===
  !  Read the slab depth from McCrory et al. (2007,2012)
  !=== 
  call read_slab_xyz()
  write(*,*) 'read slab depth done'

  !===
  !  Read the wUS-SH-2010 model 
  !===
  call read_wUS_xyz()  

  !===
  !  Generate the thermal profile at age zero
  !===
  call setupzeroage(zeroaget,zerocoeff,pottemp,adbatgrad,srate)   
  write(*,*) 'Initialize thermal model'
  !===
  !  Main part: generate model
  !=== 
  write(11,111) beglat,endlat,latinc,beglon,endlon,loninc
111   format(2(F9.3,2x,F9.3,2x,F6.3))

  !===
  ! need to find which model regime point is in.  
  ! nregion:
  ! 0 - default continental model; 
  ! 1 - oceanic model; 
  ! 2 - forearc model with slab; 
  ! 3 - normal continental model;
  !===
  write(*,*) beglat,endlat,latinc,beglon,endlon,loninc
  n_lat = nint( (endlat - beglat)/latinc ) + 1
  n_lon = nint( (endlon - beglon)/loninc ) + 1

  if(n_lat.lt.1 .or. n_lon.lt.1) then
    write(*,*) n_lat,n_lon
    write(*,*) 'check the lon/lat range!'
    stop
  endif

  do ii = 1, n_lat
    do jj = 1, n_lon
      slat = (ii-1)*latinc + beglat
      slon = (jj-1)*loninc + beglon         

      nregion = 0
      couldbe_seamnt = 0
      !write(*,*)slat,slon
      !===
      ! first, look for offshore (100 m water or deeper), from crude outline
      !===
      if ((slat.gt.coastlat(ncoast)).or.(slat.lt.coastlat(1))) then
        write(*,'(F9.3,F9.3,A40)') slat, slon, 'location out of research range'
        cycle
      endif

      !===
      ! test for location between profiles and then for location west of coast
      !===
      do ic = 1, ncoast-1
        if ((slat.ge.coastlat(ic)).and.(slat.lt.coastlat(ic+1))) then 
         
          !=== region 1: oceanic model with marine sediments ===
          ! Note: 
          !   if water is too shallow or crustal age is not available 
          !   move this point to region 2
          !===
          coint = (slat-coastlat(ic))/(coastlat(ic+1)-coastlat(ic))
          colonint = coint*coastlon(ic+1) + (1.0-coint)*coastlon(ic)

          if (slon.le.colonint) then
            if(slon.le.colonint-3) then 
              couldbe_seamnt = 1
            endif
            nregion = 1  
            call oceancrustgen(slat,slon,thickl,density,vp,vs,zmoho,&
                nmantle,sed_vs)
            if(nregion == 1) then
                write(18,*) slon,slat,nregion
                write(24,*) slon,slat,sed_vs
            endif
            exit
          endif

        endif
      enddo

      !=== test ===
      !cycle
      !=== end ===

      !===
      ! If location not off the coast, estimate moho thickness from receiver 
      ! function values.   In this version, use weighted average of points 
      ! within 0.75 or if none within 0.75 deg, within 2.0 deg weights range 
      ! from 1.0 at site to 0 at limit distance, inversely proportional to 
      ! distance
      !===

      if (nregion == 0 .or. nregion == 2) then

        call read_slab(slon,slat,slabdepth)

        !=== region 2: forearc model with slab ===
        ! Note: Slab deeper than slabmohomax (38 km) will be assign to region 3
        !===

        write(*,*) 'slabdepth',slabdepth

        if (slabdepth.lt.slabmohomax .and. slabdepth.gt.0.0) then
          nregion = 2
          write(*,*) 'forearc crust',slat,slon
          call forearccrustgen(slat,slon,thickl,density,vp,vs,zmoho,nmantle,&
                defthick,defdens,defvp,defvs,ndef)
          write(18,*)slat,slon,nregion
        else
          !=== continental model ===
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
            nregion = 3
            write(18,*)slat,slon,nregion
          else
            do i = 1,nmoho
              dist2 = (slat-moholat(i))**2 + (slon-moholon(i))**2
              dist = sqrt(dist2)
              if (dist.lt.2.0) then
                wgti = 2.0-dist
                wgtsum = wgtsum + wgti
                thicksum = thicksum + mdepth(i)*wgti
              endif
            enddo

            if (wgtsum.ne.0.0) then
              zmoho = thicksum/wgtsum
              nregion = 3
              write(18,*)slat,slon,nregion
            else if(slon >= -127.0) then
              zmoho = slabmohomax
              nregion = 3
              write(18,*)slat,slon,nregion
            else
              ! If too far from any receiver-function-constrained point, 
              ! use 31 km crust.
              zmoho = 31.0
              nregion = 0
              write(18,*)slat,slon,nregion
            endif
          endif

        endif
          
        !===
        ! If point not within any other region, use default crustal model 
        ! and receiver function crustal thickness
        !===
        if (nregion.eq.0 .or. nregion.eq.3) then
          call defaultgen(slat,slon,zmoho,defthick,defdens,defvp,defvs,ndef)
        endif

      endif

      write(17,300) slon,slat,zmoho
300   format (2F9.3,1x,f4.1)
    enddo
  enddo

  !===
  ! replace default mantle with Brandon's model
  !===
  close(11)
  close(23)
  close(24)
  close(19)

  open(21, file = outputfile)
  open(91, file = 'tempoutput')

  call mantlevar

  close(unit=91)
  close(unit=17)
  close(unit=18)
  close(unit=21)

  stop
end


!=== 
! Read body wave model wUS-SH-2010_percent.nc (Schmandt and Humphreys, 2010)
! Input: file 11, 22, 12
! Output: file 21
!===
subroutine mantlevar
  use wUSpara
  real :: defthick(50), defdens(50), defvp(50), defvs(50), dpth(50)
  real :: defthick2(50),defdens2(50),defvp2(50),defvs2(50)
  real :: blat(100),blon(100),bdepth(30),bdVs(100,100,30)
  real :: minlat,minlon,maxlat,maxlon,latinc,mindepth,maxdepth

  real, parameter :: vs2vp = 0.7
  real, parameter :: vs2dens = 0.0

  integer :: nx, ny, nxy
  integer :: i, j, k, ii, inxy
  
  read(91,*) beglat,endlat,dlat,beglon,endlon,dlon
  !read(22,*) beglat,endlat,dlat, beglon,endlon,dlon
  write(21,120) beglat,endlat,dlat,beglon,endlon,dlon
120   format(2(F9.3,1x,F9.3,1x,F6.2))

  !===
  !  11 is ab initio starting model
  !  22 is starting model using result of previous inversion.  Need both
  !  because Brandon's model is percent variation from a uniform velocity
  !  mantle model, so need uniform mantle from 11, but need new crustal 
  !  results from 22
  !===

  nx = nint( (endlon-beglon)/dlon + 1)
  ny = nint( (endlat-beglat)/dlat + 1)
  nxy =  nx * ny

  !read(13,*) minlat,minlon,maxlat,maxlon,latinc,mindepth,&
  !           maxdepth,ndepth,vs2vp,vs2dens
  !===
  ! these parameters describe range of body wave grid.  vs2vp and vs2dens describe
  ! ratio of p and density perturbations to shear wave perturbations (%) - 
  ! vs3vp = 1 keeps constant poisson's ratio
  !===

  !nby = (maxlat-minlat)/latinc + 1.01
  !nbx = (maxlon-minlon)/latinc + 1.01
  !nbxy = nbx*nby

  !===
  !  this routine for introducing perturbations to uniform mantle model assumes
  !  that grid spacing same for mantle model and perturbation model. although 
  !  there may be fewer points in perturbation model.  Also, same number of depths
  !  for every lat,long point in Brandon's model.
  !===
  !do j = 1, nbx
  !  do i = 1,nby
  !    do k = 1, ndepth
  !      read(13,*) blat(i),blon(j),bdepth(k),bdVs(i,j,k)
  !    enddo
  !  enddo
  !enddo
        

  do inxy = 1, nxy
    read(91,*,iostat=io) slon,slat
    if (io/=0) exit
    read(91,*) nlay

    write(21,*) slon,slat
    write(21,*) nlay     

    nnx = nint( (slon-beglon)/dlon + 1)
    nny = nint( (slat-beglat)/dlat + 1)
    nnidx = (nny - 1)*nx + nnx
 
    cumlayer = 0.0
    do ii = 1, nlay
      read(91,*,iostat=io) defthick(ii),defdens(ii),defvp(ii),defvs(ii)
      if(defthick(ii).lt.0.4) write(*,*) slon,slat,ii,defthick(ii),&
          'layer too thin!'
      !read(22,*) defthick2(ii),defdens2(ii),defvp2(ii),defvs2(ii)
      dpth(ii) = cumlayer + 0.5*defthick(ii)
      cumlayer = cumlayer + defthick(ii)
    enddo

    !write(*,*)slon,slat,nlay,nnx,nny,nnidx

    !=== check whether within body wave grid ===
    if  (((slon.ge.wus_lon_min).and.(slon.le.wus_lon_max)).and. &
         ((slat.ge.wus_lat_min).and.(slat.le.wus_lat_max))) then

      i = nint((slat-wus_lat_min)/wus_lat_inc + 1)
      j = nint((slon-wus_lon_min)/wus_lat_inc + 1)

      do ii = 1,nlay

        ! check whether in crust or shallower (or deeper) than depth 
        ! limits of body wave model
        if ((dpth(ii).lt.wus_dep_min).or.(dpth(ii).gt.wus_dep_max).or.&
            &(defdens(ii).lt.3.2)) then
          write(21,119) defthick(ii),defdens(ii),defvp(ii),defvs(ii)
119   format(F7.3,2x, 3(F6.3,3x))
          !write(21,*) defthick2(ii),defdens2(ii),defvp2(ii),defvs2(ii)
        else

          ! Find depth interval, then linearly interpolate perturbation 
          ! value to middle of layer. If top or bottom, interpolate to 
          ! zero at mindepth or maxdepth
           if (dpth(ii).lt.wusdep(1)) then
             prtrbfrac = 0.01*wusdvs(i,j,1)*(dpth(ii)-wus_dep_min)/&
                 (wusdep(1)-wus_dep_min)
           endif

           if (dpth(ii).gt.wusdep(ndep)) then
             prtrbfrac = 0.01*wusdvs(i,j,ndep)*(wus_dep_max-dpth(ii))/&
                 (wus_dep_max-wusdep(ndep))
           endif

           if ((dpth(ii).ge.wusdep(1)).and.(dpth(ii).le.wusdep(ndep))) then
             do jj = 1, ndep-1
               if ((dpth(ii).ge.wusdep(jj)).and.(dpth(ii).lt.wusdep(jj+1)))then 
                 inddpth = jj
                 exit
               endif
             enddo

             prtrbfrac=0.01*(wusdvs(i,j,inddpth) + (wusdvs(i,j,inddpth+1) - &
                 wusdvs(i,j,inddpth))* (dpth(ii)-wusdep(inddpth)) /&
                 (wusdep(inddpth+1)-wusdep(inddpth)))
           endif

           defvs(ii)  = (1.0+prtrbfrac)*defvs(ii)
           defvp(ii)  = (1.0+vs2vp*prtrbfrac)*defvp(ii)
           defdens(ii)= (1.0+vs2dens*prtrbfrac)*defdens(ii)
!          write(*,*) slat,slon,i,j,prtrbfrac
!          write(*,*) dpth(ii),defthick(ii),defdens(ii),defvp(ii),defvs(ii)
           write(21,119) defthick(ii),defdens(ii),defvp(ii),defvs(ii)
         endif
       enddo
    else
      do ii = 1, nlay
        !write(21,*) defthick2(ii),defdens2(ii),defvp2(ii),defvs2(ii)
        write(21,119) defthick(ii),defdens(ii),defvp(ii),defvs(ii)
      enddo
    endif

  enddo

  return
  end

  
!====================================================================
! Generate default continetal crust and mantle model 
! read moho depth and adjust layers in crust and mantle accordingly
!====================================================================
subroutine defaultgen(slat,slon,zmoho,defthick,dfdens,dfvp,dfvs,ndef)
  use region
  use topopara
  real :: defthick(50),defdens(50),defvp(50),defvs(50),dpth(50)
  real :: dfdens(50),dfvs(50),dfvp(50)
  real :: drho(3),dvp(3),dvs(3)
  real :: thick_merge,dens_merge,vp_merge,vs_merge
  real :: mthick(50),mdens(50),mvp(50),mvs(50)

  real :: thick1, thick2
  real :: slon, slat, zmoho
  integer :: i_next
  integer :: n_model
  integer :: waterd
  real :: waterdp
  ! find layer that moho lies in - make 1 km the minimum layer thickness
  ! This routine specific to default crust with 6 layers

  do i = 1,ndef
    defdens(i) = dfdens(i)
    defvs(i)   = dfvs(i)
    defvp(i)   = dfvp(i)
  enddo

  n_model = 0

  !===
  ! modify crustal properties upper 3 layers based on 8-sec Rayleigh wave
  ! dispersion data from Moschetti et al. (2010)
  !===
  call moschetti(slat,slon,drho,dvp,dvs)

  do i = 1,3
    defdens(i) = defdens(i) + drho(i)
    defvs(i)   = defvs(i)   + dvs(i)
    defvp(i)   = defvp(i)   + dvp(i)
  enddo  
    
  depth = 0.0

  do i = 1, ndef
    depth   = depth + defthick(i)
    dpth(i) = depth
  enddo

  do i = 2,ndef
    if ((zmoho.gt.dpth(i-1)).and.(zmoho.le.dpth(i))) then
       idepth = i
       !write(*,*)'Moho at depth',zmoho,'in layer',idepth
    endif
  enddo

  !===
  !  if moho lies within original 6 layers, truncate crust and extend first
  !  mantle layer properties upward , i.e. layer 7)
  !===
  if (idepth.le.6) then

    thick1 = zmoho - dpth(idepth-1)
    thick2 = dpth(idepth) - zmoho

    do i = 1, idepth-2
      !write(11,100) defthick(i),defdens(i),defvp(i),defvs(i)
      n_model = n_model + 1
      mthick(n_model) = defthick(i)
      mdens(n_model)  = defdens(i)
      mvp(n_model)    = defvp(i)
      mvs(n_model)    = defvs(i)  
    enddo

    !=== insert moho ===
    ! 1, merge thin layer 
    if (thick1.lt.1.0) then
      call mergelay(defthick(idepth-1), defdens(idepth-1), defvp(idepth-1),&
          defvs(idepth-1), thick1, defdens(idepth), defvp(idepth), &
          defvs(idepth), thick_merge, dens_merge, vp_merge, vs_merge)

      !write(11,100) thick_merge,dens_merge,vp_merge,vs_merge
      !write(11,100) thick2,defdens(7),defvp(7),defvs(7)
      n_model = n_model + 1
      mthick(n_model) = thick_merge
      mdens(n_model)  = dens_merge
      mvp(n_model)    = vp_merge
      mvs(n_model)    = vs_merge

      n_model = n_model + 1
      mthick(n_model) = thick2
      mdens(n_model)  = defdens(7)
      mvp(n_model)    = defvp(7)
      mvs(n_model)    = defvs(7)
      
      i_next = idepth + 1
    ! 2, always fill the space under moho with mantle structure
    else if (thick2.lt.1.0) then
      thick2  = thick2 + defthick(idepth+1)
      !write(11,100) defthick(idepth-1),defdens(idepth-1),defvp(idepth-1),defvs(idepth-1)
      n_model = n_model + 1
      mthick(n_model) = defthick(idepth-1)
      mdens(n_model)  = defdens(idepth-1)
      mvp(n_model)    = defvp(idepth-1)
      mvs(n_model)    = defvs(idepth-1)
      !write(11,100) thick1,defdens(idepth),defvp(idepth),defvs(idepth)
      n_model = n_model + 1
      mthick(n_model) = thick1
      mdens(n_model)  = defdens(idepth)
      mvp(n_model)    = defvp(idepth)
      mvs(n_model)    = defvs(idepth)
      !write(11,100) thick2,defdens(7),defvp(7),defvs(7)
      n_model = n_model + 1
      mthick(n_model) = thick2
      mdens(n_model)  = defdens(7)
      mvp(n_model)    = defvp(7)
      mvs(n_model)    = defvs(7)

      i_next = idepth + 2
    else
      !write(11,100) defthick(idepth-1),defdens(idepth-1),defvp(idepth-1),defvs(idepth-1)
      n_model = n_model + 1
      mthick(n_model) = defthick(idepth-1)
      mdens(n_model)  = defdens(idepth-1)
      mvp(n_model)    = defvp(idepth-1)
      mvs(n_model)    = defvs(idepth-1)
      !write(11,100) thick1,defdens(idepth),defvp(idepth),defvs(idepth)
      n_model = n_model + 1
      mthick(n_model) = thick1
      mdens(n_model)  = defdens(idepth)
      mvp(n_model)    = defvp(idepth)
      mvs(n_model)    = defvs(idepth)
      !write(11,100) thick2,defdens(7),defvp(7),defvs(7)
      n_model = n_model + 1
      mthick(n_model) = thick2
      mdens(n_model)  = defdens(7)
      mvp(n_model)    = defvp(7)
      mvs(n_model)    = defvs(7)

      i_next = idepth + 1
    endif
    

    if (i_next.le.6) then
      do i = i_next, 6
        !write(11,100) defthick(i),defdens(7),defvp(7),defvs(7)
        n_model = n_model + 1
        mthick(n_model) = defthick(i)
        mdens(n_model)  = defdens(7)
        mvp(n_model)    = defvp(7)
        mvs(n_model)    = defvs(7)
      enddo

      i_next = 7
    endif

  endif

  !=== 
  ! if moho lies in the default mantle, extend the bottom layer of crust. 
  !===
  if (idepth.gt.6) then

    thick1 = zmoho - dpth(idepth-1)
    thick2 = dpth(idepth) - zmoho

    ! 1, copy the top crustal layers
    do i = 1, 5
      !write(11,100) defthick(i),defdens(i),defvp(i),defvs(i)
      n_model = n_model + 1
      mthick(n_model) = defthick(i)
      mdens(n_model)  = defdens(i)
      mvp(n_model)    = defvp(i)
      mvs(n_model)    = defvs(i)  
    enddo

    ! 2, extend bottom layer when necessary 
    if (idepth.gt.7) then
      do i = 6,idepth-2
        !write(11,100) defthick(i),defdens(6),defvp(6),defvs(6)
        n_model = n_model + 1
        mthick(n_model) = defthick(i)
        mdens(n_model)  = defdens(6)
        mvp(n_model)    = defvp(6)
        mvs(n_model)    = defvs(6)  
      enddo
    endif

    !=== insert moho ===
    !1, adjust thin layer 
    if (thick1.lt.2.0) then
      thick_merge = thick1 + defthick(idepth-1)
      !write(11,100) thick1,defdens(6),defvp(6),defvs(6)
      n_model = n_model + 1
      mthick(n_model) = thick_merge
      mdens(n_model)  = defdens(6)
      mvp(n_model)    = defvp(6)
      mvs(n_model)    = defvs(6)
      !write(11,100) thick2,defdens(idepth),defvp(idepth),defvs(idepth)
      n_model = n_model + 1
      mthick(n_model) = thick2
      mdens(n_model)  = defdens(idepth)
      mvp(n_model)    = defvp(idepth)
      mvs(n_model)    = defvs(idepth)

      i_next = idepth + 1 
    else if (thick2.lt.2.0) then
      call mergelay(defthick(idepth+1),defdens(idepth+1),defvp(idepth+1),defvs(idepth+1),&
                    thick2,            defdens(idepth),  defvp(idepth),  defvs(idepth),&
                    thick_merge,          dens_merge,       vp_merge,       vs_merge)

      !write(11,100) defthick(idepth-1),defdens(6),defvp(6),defvs(6)
      n_model = n_model + 1
      mthick(n_model) = defthick(idepth-1)
      mdens(n_model)  = defdens(6)
      mvp(n_model)    = defvp(6)
      mvs(n_model)    = defvs(6)
      !write(11,100) thick1,defdens(6),defvp(6),defvs(6)
      n_model = n_model + 1
      mthick(n_model) = thick1
      mdens(n_model)  = defdens(6)
      mvp(n_model)    = defvp(6)
      mvs(n_model)    = defvs(6)
      !write(11,100) thick2,dens_merge,vp_merge,vs_merge
      n_model = n_model + 1
      mthick(n_model) = thick_merge
      mdens(n_model)  = dens_merge
      mvp(n_model)    = vp_merge
      mvs(n_model)    = vs_merge

      i_next = idepth + 2
    else 
      !write(11,100) defthick(idepth-1),defdens(6),defvp(6),defvs(6)
      n_model = n_model + 1
      mthick(n_model) = defthick(idepth-1)
      mdens(n_model)  = defdens(6)
      mvp(n_model)    = defvp(6)
      mvs(n_model)    = defvs(6)
      !write(11,100) thick1,defdens(6),defvp(6),defvs(6)
      n_model = n_model + 1
      mthick(n_model) = thick1
      mdens(n_model)  = defdens(6)
      mvp(n_model)    = defvp(6)
      mvs(n_model)    = defvs(6)
      !write(11,100) thick2,defdens(idepth),defvp(idepth),defvs(idepth)
      n_model = n_model + 1
      mthick(n_model) = thick2
      mdens(n_model)  = defdens(idepth)
      mvp(n_model)    = defvp(idepth)
      mvs(n_model)    = defvs(idepth)

      i_next = idepth + 1
    endif

  endif
   
  do i = i_next, ndef
    !write(11,100) defthick(i),defdens(i),defvp(i),defvs(i)
    n_model = n_model + 1
    mthick(n_model) = defthick(i)
    mdens(n_model)  = defdens(i)
    mvp(n_model)    = defvp(i)
    mvs(n_model)    = defvs(i)  
  enddo

  write(11,101) slon,slat
  write(11,*) n_model

  do i = 1, n_model
    write(11,100) mthick(i),mdens(i),mvp(i),mvs(i)
  enddo
100   format(F7.3,2x, 3(F6.3,3x))
101   format(F9.3,2x,F7.3)

  call readwatdep2(slon,slat,waterd)
  if (waterd < 0) then
    waterdp = waterd / 1000.
  else
    waterdp = 0.0
  endif

  write(23,301) slon,slat,waterdp
301   format (2F9.3,1x,f4.1)
  return
end        

!================================================================================
!  find points within 0.6 deg of slat,slon and average them.  Moschetti on 0.5 deg
!  grid.  Averaging provides both interpolation and smoothing - seems to be some
!  oscillation in their solution with largest values adjacent to smallest, so do a
!  little smoothing.
! 
!  background value for phase velocity in their model is 3.0772 km/s at 8 s period
!================================================================================= 
subroutine moschetti(slat,slon,drho,dvp,dvs)
  use moschettil
  use region
  real :: drho(3),dvp(3),dvs(3)
  logical :: is_mos

  weight = 0.0
  sum = 0.0
  is_mos = .false.

  do i = 1,nmos
    dist = sqrt((slat-moslat(i))**2 + (slon-moslon(i))**2)
    if (dist.le.0.6) then
      is_mos = .true.
      wgt = 1.0/(dist+0.3) 
      weight = weight + wgt
      sum = sum + mosrayl8(i)*wgt
    endif
  enddo

  if(is_mos) then
    davg = sum/weight - 3.0772
  else 
    davg = 0.0
  endif

  !  limit variability of model
  if (davg.gt.0.2)    davg = 0.2
  if (davg.lt.(-0.4)) davg = -0.4

  !  get about 0.1 km/s change in phase velocity for 0.4 km/s change in Vs in upper
  !  3 km, 0.2 km/s in next 3 km and 0.1 in next 4.  Not really linear, but assume
  !  it is, with the limitations of linearity assumption ameliorated by the
  !  limitation on variability in the model.  (This assumes dVp = dVs and
  !  drho = dVp/2.2

   dvs(1) = 4.0*davg
   dvs(2) = 2.0*davg
   dvs(3) = 1.0*davg

   do i = 1,3
     dvp(i) = dvs(i)
     drho(i)= dvs(i)/2.2
   enddo

   return
end

!===============================================================================
!  generate crustal layers and mantle layers down to 210 km - crustal layering
!  based on Lin,Shearer etal. tomo and Levanders's moho - usedfixed lower crust 
!  velocity below 22 km, since is not constrained anyway by tomo
!
!  crustal layers are 0-3, 3-6, 6-10, 10-15, 15-22, 22-moho
!  assumes mantle layers 10 km thick down to 90 km
!===============================================================================
subroutine socalgen(slat,slon,zmoho,imin,jmin,slvp,slvpvs,thickl,density,vp,vs,nmantle)
  use region
  real :: slvp(29,45,8),slvpvs(29,45,8)
  real :: thick(5),cvp(5),cvs(5),cden(5)
  real :: thickl(50),density(50),vp(50),vs(50)


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
 
    if (zmoho.gt.50.0) write(*,'(3F9.3,A20)') slat,slon, zmoho,' crust over 50 km'

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

  if (zmoho.lt.16.0) write(*,'(3F9.3,A20)') slat,slon,zmoho,' crust under 16 km'

  do i = nst, nmantle
    write(11,100) thickl(i),density(i),vp(i),vs(i)
  enddo
100   format(F7.3,2x, 3(F6.3,3x))
101   format(F9.3,2x,F7.3)
  return
end
    
     
!=============================================================================  
!  Notes from D. W. Forsyth:
!  generate crustal layers and mantle layers down to 210 km - 
!  crustal layering based on waterdepth alone.  
!  Water layer directly interpolated value 
!  If water depth >= than 3000 m, then use standard oceanic crust with 
!  2 km upper crust P velocity 5.0 km/s  
!  5 km lower crust P velocity 6.8 km/s
!
!  Interpolate shallower depths to at 100m water, 
!  constant uppermost crust, 
!  6.0 km/s layer varying from 0 to 13 km thickness, and 
!  6.8 km/s layer varying from 5 to 13 km thickness.
!
!  assumes mantle layers 10 km thick down to 90 km (beginning at 10 km)
!  if crust thickens, thin or eliminate mantle layers to compensate 
!  
!  === YYR 11/19/14 Princeton ===
!  For JdF case:
!   sediment when thicker than 0.2 km will be included
!   crust 2A layer with constant thickness 0.4 km is included
!       
!  
!==============================================================================
subroutine oceancrustgen(slat,slon,thickl,density,vp,vs,zmoho,nmantle,sed_vs)
  use globalage
  use region
  use topopara
  use thermodel
  integer, parameter :: nmaxly = 50
  real :: thickl(nmaxly),density(nmaxly),vp(nmaxly),vs(nmaxly)

  real, parameter :: waterden = 1.03
  real, parameter :: watersv  = 0.00
  real, parameter :: waterpv  = 1.50

  real, parameter :: rhosed   = 2.00

  real, parameter :: thick2A  = 0.40

  !=== default layer 2B and 3 in crust ===
  !real, parameter :: thickup  = 1.60
  real, parameter :: thickup  = 1.60
  real, parameter :: denup    = 2.45
  real, parameter :: svup     = 2.63
  real, parameter :: pvup     = 5.00

  real, parameter :: thicklow = 5.00
  real, parameter :: denlow   = 3.05
  real, parameter :: svlow    = 3.89
  real, parameter :: pvlow    = 6.80

  !=== sea floor age ===
  real, parameter :: cerr = -999.0
  real, parameter :: epsiln = 1.0E-06
  real, parameter :: m_age = 8.0 ! constant age for mantle structure
  real :: age, diffage

  !=== water depth and sediment thickness ===
  integer :: waterd, sedthck
  real :: waterdp, hsed, halfhsed, sed_vs

  !=== layer 2A ===
  real :: vp2A, rho2A, vs2A

  logical :: is_sed, is_wat, is_age

  real :: thick_merge,dens_merge,vp_merge,vs_merge
  real :: m_depth, vs_t, vp_t

  integer :: i 
  !=== water layer ===
  ! read topo data from ETOPO_1
  call readwatdep2(slon,slat,waterd)
  !write(*,*)'water depth:', waterd
  !if(waterd < 100. .and. waterd > 0.) then
  if(waterd < 400. .and. waterd > 0.) then
    if(couldbe_seamnt == 1) then
      is_wat = .true.  
      waterd = 400
      waterdp = waterd/1000.  
    else
      is_wat = .false.
      waterdp = 0.0
      nregion = 2 ! move to forearc model
      write(*,'(2F9.3,A9,I6,A30,I2)') slon,slat,'water =',waterd,'move to forearc: 1 ->',nregion
      return
     endif  
  else if (waterd < 0.) then
    ! seamount ? 
    is_wat = .false.
    nregion = 0 ! move to continental model
    write(*,'(2F9.3,A10,I6,A30,I2)') slon,slat,'water =',waterd,'move to continent: 1 ->',nregion
    return
  else
    is_wat = .true.
    waterdp = waterd/1000. ! in km
  endif


  !=== marine sediment layer ===
  ! Vp and density from Hamilton (1979)
  ! Vs from Ruan et al. (2014)
  !======================
  call read_sed(slon,slat,sedthck)
  !write(*,*)'sediment thickness:', sedthck

  ! add shallow water to marine sediments
  !if(is_wat == .false. .and. waterd > 0) then
  !  sedthck = sedthck + waterd
  !endif 

  ! if marine sediments still too thin, neglect it.
  !if(sedthck < 100.) then 
  if(sedthck < 0.) then 
    is_sed = .false.
    hsed = 0.0
    write(*,*)'sedthck',sedthck
  else
    hsed = sedthck/1000.  ! in km
    is_sed = .true.
    halfhsed = hsed/2.0   
    vssed = (sedthck + 139.3560)/1373.9115 ! fit from Ruan et al. (2014)
    vpsed = 1.511 + 1.304*halfhsed - 0.741*halfhsed**2 + 0.257*halfhsed**3
  endif

  !=== output for plotting
  sed_vs = vssed

  !=== crust layer 2A ===
  is_age = .true.

  call crustage(slon,slat,age)

  diffage = abs(age - cerr)
  if (diffage .le. epsiln) then
     call crustageglb(slon, slat, age)
  endif 
  !write(*,*)'crustal age:',age
  
  if(age .ge. 0.0) then
     call crust2A_age2vel(age,vp2A,vs2A,rho2A)
  else 
     is_age = .false.
     write(*,'(3F9.3,A30)') slon,slat,age,' no age data replace 2A with 2B'
     !nregion = 2
     !return
  endif

  !write(*,*)'Vp in layer 2A:',vp2A
  !write(*,*)'rho in layer 2A:', rho2A

100   format(F7.3,2x, 3(F6.3,3x))
101   format(F9.3,2x,F7.3)

  !if(age .lt. 0.0) then
  !  write(*,*) slon,slat,waterd,sedthck,age
  !endif

   zmoho = waterdp+hsed+thick2A+thickup+thicklow

   ! since 2A and sediments will be merged together
   if(is_sed) then
      nlayer  = nmantle + 4
   else
      nlayer  = nmantle + 3
   endif

   !if(is_age == .false.) then
   !  zmoho = waterdp+hsed+thickup+thicklow
   !  nlayer = nlayer - 1 
   !endif

   thickm1 = thickl(1) - zmoho + 10.0
  
   ! === build the model 
   write(11,101) slon, slat
   write(11,*) nlayer
   write(11,100) waterdp,waterden,waterpv,watersv

   if(is_sed) then
      if(is_age) then        
        call mergelay(hsed,rhosed,vpsed,vssed,&
                      thick2A,rho2A,vp2A,vs2A,&
                      thick_merge,dens_merge,vp_merge,vs_merge)
      else
        call mergelay(hsed,rhosed,vpsed,vssed,&
                      thick2A,denup,pvup,svup,&
                      thick_merge,dens_merge,vp_merge,vs_merge)
        
      endif  
      write(11,100) thick_merge,dens_merge,vp_merge,vs_merge
   else
      if(is_age) then
        write(11,100) thick2A,rho2A,vp2A,vs2A
      else
        write(11,100) thick2A,denup,pvup,svup  
      endif
   endif

   write(11,100) thickup,denup,pvup,svup
   write(11,100) thicklow,denlow,pvlow,svlow

   m_depth = zmoho
   call thermantle(age,m_depth,thickm1,vs_t,vp_t)
   !write(*,*) age,m_depth,vs_t,vp_t
 
   write(11,100) thickm1,density(1),vp_t,vs_t
   m_depth = m_depth + thickm1
  
   nst = 2
   do i = nst, nmantle
     if(m_depth.lt.180.and.is_age) then ! max thermal model is 200 km when age is available
       call thermantle(age,m_depth,thickl(i),vs_t,vp_t)
       !write(*,*) age,m_depth,vs_t,vp_t
       write(11,100) thickl(i),density(i),vp_t,vs_t
       m_depth = m_depth + thickl(i)
     else
       write(11,100) thickl(i),density(i),vp(i),vs(i)
     endif
   enddo

   write(23,302) slon,slat,waterdp ! liquid-solid interface
!   write(19,302) slon,slat,waterdp+hsed
!   Y Ruan 02/16/2016. comment out to generate slab surface only at continental
!   side and output water surface instead for SPECFEM3D meshing
!   write(19,302) slon,slat,waterdp ! top of slab
   write(25,302) slon, slat, 0.0
302   format (2F9.3,1x,f4.1) 
   return
end


!=============================================================================  
!  generate crustal layers and mantle layers down to 210 km - 
!
!  Interpolate shallower depths to at 100m water,
!
!  Fill in default continental crust model between seafloor (if there is water) 
!  and slab (if there is slab) 
!
!  Use standard oceanic crust with 
!  2 km upper crust P velocity 5.0 km/s  
!  5 km lower crust P velocity 6.8 km/s
!
!  use default mantle at 40 to 60 km 
!  assumes mantle layers 10 km thick down to 90 km (beginning at 10 km)
!  if crust thickens, thin or eliminate mantle layers to compensate 
!  
!  Below 60 km use schmandt and humphreys model
!
!  YYR 11/19/14 Princeton 
!===============================================================================


subroutine forearccrustgen(slat,slon,&
                           thickl,density,vp,vs,zmoho,nmantle,&
                           defthick,defdens,defvp,defvs,ndef)
  use region
  use thermodel
  use slabpara
  use topopara
  real :: slat,slon
  integer, parameter :: nmaxly = 50
  integer :: ndef
  real :: thickl(nmaxly),density(nmaxly),vp(nmaxly),vs(nmaxly)
  real :: defthick(50),defdens(50),defvp(50),defvs(50),dpth(50)

  !=== default water ===
  real, parameter :: waterden = 1.03
  real, parameter :: watersv  = 0.00
  real, parameter :: waterpv  = 1.50

  !=== default sediment ===
  real, parameter :: rhosed   = 2.00

  !=== default oceanic crust ===
  ! replace 2A with 2 to remove the unrealistic low vel layer
  real, parameter :: thickup  = 2.00
  !real, parameter :: denup    = 2.45
  !real, parameter :: svup     = 2.63
  !real, parameter :: pvup     = 5.00
  real, parameter :: denup    = 3.05
  real, parameter :: svup     = 3.89
  real, parameter :: pvup     = 6.80

  real, parameter :: thicklow = 5.00
  real, parameter :: denlow   = 3.05
  real, parameter :: svlow    = 3.89
  real, parameter :: pvlow    = 6.80

  real, parameter :: forearc_age = 8.0 ! constant age for thermal mantle
  !=== water depth and sediment thickness ===
  integer :: waterd, sedthck
  real :: waterdp, hsed, halfhsed
 
  !=== slab depth ===
  real :: slabdepth, diffslab, basement

  real :: zmoho
  real :: m_depth(10), thickm_top
  integer :: imantle_start

  !=== for added continental crust ===
  integer :: nadd, nmergetop, nmergebot
  integer :: i, iadd_start, iadd_end
  real    :: thickadd(6),densadd(6),vpadd(6),vsadd(6)

  real :: depth, man_dep, vs_t, vp_t
  logical :: is_sed, is_wat, is_slab, is_continent

  real :: ztopo
  !===
  !  Step 1: Find basement of seafloor (includes sed)
  !===

  !=== water layer ==
  is_wat = .true.
  is_sed = .true.
  is_continent = .true.
  
  call readwatdep2(slon,slat,waterd)
  if(waterd .le. 400.) then
    is_wat = .false.
    waterdp = 0.0
  else
    waterdp = waterd/1000. ! in km
  endif

  !=== sediment layer ===
  ! Vp and density from Hamilton (1976)
  ! Vs from Ruan et al. (2014)
  !======================
  call read_sed(slon,slat,sedthck)
  if(is_wat==.false..and.waterd.gt.0.) then
    sedthick = sedthick + waterd
    write(*,*) 'replace thin water layer with sediment',waterd/1000.
  endif  

  if(sedthck < 400.) then 
    is_sed = .false.
    hsed = 0.0
  else
    hsed = sedthck/1000.  ! in km
    halfhsed = hsed/2.0   
    vssed = (sedthck + 139.3560)/1373.9115 ! fit from Ruan et al. (2014)
    vpsed = 1.511 + 1.304*halfhsed - 0.741*halfhsed**2 + 0.257*halfhsed**3
  endif

  !write(*,*)'forearc Hwat Hsed', waterdp, hsed
  basement = waterdp + hsed
  write(*,*)'forearc basement depth (km)',basement

  !===
  !  Step 2: Find the top surface of oceanic slab 
  !===

  !=== 
  ! read slab depth and figure out how to build wedge with default mantle model
  !===
  call read_slab(slon,slat,slabdepth)

  write(*,*)'forearc slab top depth (km)',slabdepth

  if(slabdepth < 0.0 .or. slabdepth > 120.) then ! no slab
    slabdepth = basement
    is_continent = .false.
  else 
    diffslab = slabdepth - basement
    !=== 
    ! Fill the continental shelf wedge with sediments or water if
    ! difference between slabdepth and (water + sed) is small
    !=== 
    if (diffslab < 0) then ! just in case
      is_continent = .false.
      slabdepth = basement
    else if(diffslab > 0 .and. diffslab < 0.1) then ! shallow slab
      is_continent = .false.
      ! chukren, check here probably add diffslab twice
      if(is_sed) then
        hsed = hsed + diffslab
      else if (is_wat) then
        waterdp = waterdp + diffslab
      else
        write(*,*) 'something wrong! region:',nregion
        return
      endif
      basement = waterdp + hsed
      slabdepth = basement
    else ! deep slab, diffslab > 0.1
      is_continent = .true.
    endif
  endif

  !===
  ! Note:
  ! is_continent = .false. => default oceanic model
  ! is_continent = .true.  => oceanic model + continental
  !===

100   format(F7.3,2x, 3(F6.3,3x))
101   format(F9.3,2x,F7.3)
  
  zmoho = slabdepth + thickup + thicklow

  write(*,*)'forearc   moho   depth (km)',zmoho
  !write(*,'(I2,A10,F7.3,A12,F7.3,A8,F7.3)') nregion,'basement=',basement,'slabdepth=',slabdepth,'zmoho=',zmoho

  !===
  !   Step 3: determin how many layers we need for model 
  !===

  !=== 
  !  See how many mantle layers we need to add
  !===

  depth = 10.0
  do i = 1, 6
    depth = depth + thickl(i)
    m_depth(i) = depth
  enddo  

  do i = 1, 6
    if(zmoho.le.m_depth(1)) then
      mdepth = 1
      exit
    endif

    if ( i > 1 .and. (zmoho.gt.m_depth(i-1)) .and. (zmoho.le.m_depth(i))) then
      mdepth = i
      exit
    endif
  enddo

  thickm_top = m_depth(mdepth) - zmoho
  !write(*,*)'thickness of top mantle layer',thickm_top

  imantle_start = mdepth 

  ! merge any mantle layer less than 1 km thick to the layer beneath
  if(thickm_top .lt. 1.0) then
    imantle_start = imantle_start + 1
    thickm_top = thickl(imantle_start) + thickm_top
  endif
 
  !===
  !  See if we need add continental crust and how many layers needed
  !===
  if(is_continent) then
    depth = 0.0

    do i = 1, ndef
      depth   = depth + defthick(i)
      dpth(i) = depth
    enddo

    ! First find where the slab and basement are located
    ! chukren check this case!!!!
    sdepth = 1
    bdepth = 1
    do i = 2,ndef
      if((slabdepth.gt.dpth(i-1)) .and. (slabdepth.le.dpth(i))) then
         sdepth = i
      endif

      if((basement.gt.dpth(i-1)) .and. (basement.le.dpth(i))) then
         bdepth = i
      endif
    enddo
    
    !write(*,*)'basement at layer',bdepth
    !write(*,*)'slab at layer',sdepth

    ! Second: add layers regardless and then merge thiner layer with their
    ! lower neighbor.
    nadd   = sdepth - bdepth + 1 
    iadd_start = 1
    iadd_end   = nadd

    nmergetop = 0
    nmergebot = 0
    !write(*,*)'add',nadd,'layers'
 
    ! add top layer
    thickadd(1)    = dpth(bdepth) - basement
    densadd(1)     = defdens(bdepth) 
    vpadd(1)       = defvp(bdepth)
    vsadd(1)       = defvs(bdepth) 
 
    ! add bot layer
    thickadd(nadd) = slabdepth - dpth(sdepth-1)
    densadd(nadd)  = defdens(sdepth)
    vpadd(nadd)    = defvp(sdepth)
    vsadd(nadd)    = defvs(sdepth) 

    ! add layers in between
    if(nadd > 2) then
      do i = 2, nadd - 1
        thickadd(i) = defthick(bdepth + i -1 )
        densadd(i)  = defdens(bdepth + i -1 )
        vpadd(i)    = defvp(bdepth + i -1)
        vsadd(i)    = defvs(bdepth + i -1)  
      enddo
    endif

    !write(*,*)'add continental layer',nadd

    ! merge thin layers
    if(thickadd(1) < 0.5) then
      nmergetop  = 1
      iadd_start = iadd_start + 1
      thickadd(2) = thickadd(2) + thickadd(1) 
    endif

    if(nadd > 1 .and. thickadd(nadd) < 0.5) then
      nmergebot  = 1
      iadd_end   = iadd_end - 1
      thickadd(iadd_end) = thickadd(iadd_end) + thickadd(iadd_end+1) 
    endif 
 
    nadd = nadd - nmergetop - nmergebot
    !write(*,*)'merge added layers',nadd
  endif


  !=== 
  !  Step 4: Build model 
  !===
  write(*,*)'water     layer ',is_wat
  write(*,*)'sediment  layer ',is_sed 
  write(*,*)'continent crust ',is_continent
  if(is_wat) then
    if(is_sed) then
      if(is_continent) then
        !=== 
        ! is_wat + is_sed + is_continent
        ! water + sed + continent crust + oceanic crust
        !===
        nlayer = 2 + nadd + 2 + (nmantle - imantle_start + 1)
      
        write(11,101) slon, slat
        write(11,*) nlayer
        write(11,100) waterdp, waterden,waterpv,watersv
        write(11,100) hsed, rhosed, vpsed, vssed
        do i = iadd_start, iadd_end 
          write(11,100)thickadd(i), densadd(i), vpadd(i), vsadd(i) 
        enddo
        write(11,100) thickup, denup,pvup,svup
        write(11,100) thicklow, denlow,pvlow,svlow
        
        man_dep = zmoho  
        call thermantle(forearc_age,man_dep,thickm_top,vs_t,vp_t)
        !write(*,*) forearc_age,man_dep,vs_t,vp_t
        man_dep = man_dep + thickm_top
        !write(11,100) thickm_top,density(imantle_start),vp(imantle_start),vs(imantle_start)
        write(11,100) thickm_top,density(imantle_start),vp_t,vs_t
        
        do i = imantle_start+1, nmantle
          if(man_dep .lt. 180) then
            call thermantle(forearc_age,man_dep,thickl(i),vs_t,vp_t)
            !write(*,*) forearc_age,man_dep,vs_t,vp_t
            write(11,100) thickl(i),density(i),vp_t,vs_t
            man_dep = man_dep + thickl(i)
          else
            write(11,100) thickl(i),density(i),vp(i),vs(i)
          endif
        enddo
      else 
        !=== 
        ! is_wat + is_sed 
        ! water + sed + oceanic crust
        !===
        nlayer = 2 + 2 + (nmantle - imantle_start + 1)

        write(11,101) slon, slat
        write(11,*) nlayer
        write(11,100) waterdp, waterden,waterpv,watersv
        write(11,100) hsed, rhosed, vpsed, vssed
        write(11,100) thickup, denup,pvup,svup
        write(11,100) thicklow, denlow,pvlow,svlow

        man_dep = zmoho  
        call thermantle(forearc_age,man_dep,thickm_top,vs_t,vp_t)
        !write(*,*) forearc_age,man_dep,vs_t,vp_t
        man_dep = man_dep + thickm_top
        write(11,100) thickm_top,density(imantle_start),vp_t,vs_t

        do i = imantle_start+1, nmantle
          if(man_dep .lt. 180) then
            call thermantle(forearc_age,man_dep,thickl(i),vs_t,vp_t)
            !write(*,*) forearc_age,man_dep,vs_t,vp_t
            write(11,100) thickl(i),density(i),vp_t,vs_t
            man_dep = man_dep + thickl(i)
          else 
            write(11,100) thickl(i),density(i),vp(i),vs(i)
          endif
        enddo
      endif
    else
      if(is_continent) then
        !=== 
        ! is_wat + is_continent
        ! water + continent crust + oceanic crust
        !===
        nlayer = 1 + nadd + 2 + (nmantle - imantle_start + 1)
      
        write(11,101) slon, slat
        write(11,*) nlayer
        write(11,100) waterdp, waterden,waterpv,watersv
        do i = iadd_start, iadd_end 
          write(11,100) thickadd(i), densadd(i), vpadd(i), vsadd(i) 
        enddo
        write(11,100) thickup, denup,pvup,svup
        write(11,100) thicklow, denlow,pvlow,svlow

        man_dep = zmoho  
        call thermantle(forearc_age,man_dep,thickm_top,vs_t,vp_t)
        !write(*,*) forearc_age,man_dep,vs_t,vp_t
        man_dep = man_dep + thickm_top
        write(11,100) thickm_top,density(imantle_start),vp_t,vs_t

        do i = imantle_start+1, nmantle
          if(man_dep .lt. 180) then
            call thermantle(forearc_age,man_dep,thickl(i),vs_t,vp_t)
            !write(*,*) forearc_age,man_dep,vs_t,vp_t
            write(11,100) thickl(i),density(i),vp_t,vs_t
            man_dep = man_dep + thickl(i)
          else 
            write(11,100) thickl(i),density(i),vp(i),vs(i)
          endif
        enddo
      else 
        !=== 
        ! is_wat 
        ! water + oceanic crust
        !===
        nlayer = 1 + 2 + (nmantle - imantle_start + 1)

        write(11,101) slon, slat
        write(11,*) nlayer
        write(11,100) waterdp, waterden,waterpv,watersv
        write(11,100) thickup, denup,pvup,svup
        write(11,100) thicklow, denlow,pvlow,svlow

        man_dep = zmoho  
        call thermantle(forearc_age,man_dep,thickm_top,vs_t,vp_t)
        !write(*,*) forearc_age,man_dep,vs_t,vp_t
        man_dep = man_dep + thickm_top
        write(11,100) thickm_top,density(imantle_start),vp_t,vs_t

        do i = imantle_start+1, nmantle
          if(man_dep .lt. 180) then
            call thermantle(forearc_age,man_dep,thickl(i),vs_t,vp_t)
            !write(*,*) forearc_age,man_dep,vs_t,vp_t
            write(11,100) thickl(i),density(i),vp_t,vs_t
            man_dep = man_dep + thickl(i)
          else 
            write(11,100) thickl(i),density(i),vp(i),vs(i)
          endif
        enddo

      endif
    endif
  else ! no water
    if(is_sed) then
      !===
      ! is_sed + is_continent
      ! sed + continent crust + oceanic crust
      !===
      if(is_continent) then
        nlayer = 1 + nadd + 2 + (nmantle - imantle_start + 1)
        !write(*,*)nlayer,iadd_start, iadd_end,imantle_start, nmantle
        write(11,101) slon, slat
        write(11,*) nlayer
        write(11,100) hsed, rhosed, vpsed, vssed
        do i = iadd_start, iadd_end
          write(11,100) thickadd(i), densadd(i), vpadd(i), vsadd(i) 
        enddo
        write(11,100) thickup, denup,pvup,svup
        write(11,100) thicklow, denlow,pvlow,svlow

        man_dep = zmoho  
        call thermantle(forearc_age,man_dep,thickm_top,vs_t,vp_t)
        !write(*,*) forearc_age,man_dep,vs_t,vp_t
        man_dep = man_dep + thickm_top
        write(11,100) thickm_top,density(imantle_start),vp_t,vs_t

        do i = imantle_start+1, nmantle
          if(man_dep .lt. 180) then
            call thermantle(forearc_age,man_dep,thickl(i),vs_t,vp_t)
            !write(*,*) forearc_age,man_dep,vs_t,vp_t
            write(11,100) thickl(i),density(i),vp_t,vs_t
            man_dep = man_dep + thickl(i)
          else 
            write(11,100) thickl(i),density(i),vp(i),vs(i)
          endif
        enddo
      !===
      ! is_sed 
      ! sed + oceanic crust
      !===
      else 
        nlayer = 1 + 2 + (nmantle - imantle_start + 1)
        !write(*,*)nlayer,iadd_start, iadd_end,imantle_start, nmantle
        write(11,101) slon, slat
        write(11,*) nlayer
        write(11,100) hsed, rhosed, vpsed, vssed
        write(11,100) thickup, denup,pvup,svup
        write(11,100) thicklow, denlow,pvlow,svlow

        man_dep = zmoho  
        call thermantle(forearc_age,man_dep,thickm_top,vs_t,vp_t)
        !write(*,*) forearc_age,man_dep,vs_t,vp_t
        man_dep = man_dep + thickm_top
        write(11,100) thickm_top,density(imantle_start),vp_t,vs_t

        do i = imantle_start+1, nmantle
          if(man_dep .lt. 180) then
            call thermantle(forearc_age,man_dep,thickl(i),vs_t,vp_t)
            !write(*,*) forearc_age,man_dep,vs_t,vp_t
            write(11,100) thickl(i),density(i),vp_t,vs_t
            man_dep = man_dep + thickl(i)
          else 
            write(11,100) thickl(i),density(i),vp(i),vs(i)
          endif
        enddo
      endif      
      write(*,*)'sediment w/o water, something may be odd here!'
    else
      !===
      ! is_continent
      ! continent crust + oceanic crust
      !===
      if(is_continent) then
        nlayer = nadd + 2 + (nmantle - imantle_start + 1)
        !write(*,*)nlayer,iadd_start, iadd_end,imantle_start, nmantle
        write(11,101) slon, slat
        write(11,*) nlayer
        do i = iadd_start, iadd_end 
          write(11,100) thickadd(i), densadd(i), vpadd(i), vsadd(i) 
        enddo
        write(11,100) thickup, denup,pvup,svup
        write(11,100) thicklow, denlow,pvlow,svlow

        man_dep = zmoho  
        call thermantle(forearc_age,man_dep,thickm_top,vs_t,vp_t)
        !write(*,*) forearc_age,man_dep,vs_t,vp_t
        man_dep = man_dep + thickm_top
        write(11,100) thickm_top,density(imantle_start),vp_t,vs_t

        do i = imantle_start+1, nmantle
          if(man_dep .lt. 180) then
            call thermantle(forearc_age,man_dep,thickl(i),vs_t,vp_t)
            !write(*,*) forearc_age,man_dep,vs_t,vp_t
            write(11,100) thickl(i),density(i),vp_t,vs_t
            man_dep = man_dep + thickl(i)
          else 
            write(11,100) thickl(i),density(i),vp(i),vs(i)
          endif
        enddo
      !===
      !
      ! oceanic crust 
      !===
      else 
        write(*,*) slon, slat,'check this locate!!!!!!!!!!!!!!!!'
        nlayer = 2 + (nmantle - imantle_start + 1)
        !write(*,*)nlayer,iadd_start, iadd_end,imantle_start, nmantle
        write(11,101) slon, slat
        write(11,*) nlayer
        write(11,100) thickup, denup,pvup,svup
        write(11,100) thicklow, denlow,pvlow,svlow
        write(11,100) thickm_top,density(imantle_start),vp(imantle_start),vs(imantle_start)
        do i = imantle_start+1, nmantle
          write(11,100) thickl(i),density(i),vp(i),vs(i)
        enddo
      endif      
    endif
  endif

  if (waterd < 0) then
    write(23,303) slon,slat,waterd/1000.
  else
    write(23,303) slon,slat,waterdp
  endif

  write(19,303) slon,slat,slabdepth
  ! if there is water, always output water surface
  if (is_wat) then
    write(25,303) slon,slat,0.0
  endif

303   format (2F9.3,1x,f4.1)

  return
  end


  !=== 
  !  Convert crust age to Vp velocity in layer 2A
  !  Vp = a + b*exp(c*age)
  !  from Dunn and Forsyth (2007) and Nedimovic et al. (2008)  
  !
  !  density of 2A: rho = 3.5 - 3.79/vp
  !  from Calson & Raskin (1984).
  !
  !  Poisson ratio from Collier and Singh (JGR, 1988)
  !  0.4 for 2A and 0.28 for 2A/B transition at EPR 14S
  !
  !  Average Poisson ration for oceanic crust is
  !  0.30 for 2A 
  !  0.28 for 2B
  !  0.31 for 3 (lower)
  !  0.24 for uppermost mantle 
  !  from Hyndman (Tectonophysics, 59, 321-333, 1979)
  !===
  subroutine crust2A_age2vel(age,vp,vs,rho)
  implicit none
  real, parameter :: a =  4.249
  real, parameter :: b = -1.817
  real, parameter :: c = -0.388
  real, parameter :: maxage = 42.
  real, parameter :: poisson = 0.30
  real :: age, vp, vs, rho

  if(age.ge.maxage) then
    age = maxage
    write(*,'(A12,F6.2,A36)')'Warning: crustal age=',age,'Ma is reset to 42 Ma'
    ! chukren. Modified for an enlarged area
    ! set any age older than 40 Ma to 40 Ma. 
    !stop
  endif

  vp = a + b*exp(c*age)
  vs = sqrt(1.0 - 0.5/(1.0 - poisson)) * vp
  rho = 3.5 - 3.79/vp

  return
  end

  subroutine mergelay(thick_m1,dens_m1,vp_m1,vs_m1,&
                      thick_m2,dens_m2,vp_m2,vs_m2,&
                      thick_m,dens_m,vp_m,vs_m)
  implicit none
  real, intent(in) :: thick_m1,dens_m1,vp_m1,vs_m1
  real, intent(in) :: thick_m2,dens_m2,vp_m2,vs_m2
  real :: thick_m,dens_m,vp_m,vs_m
  real :: ttime_sum,thick_sum

  thick_sum = thick_m1 + thick_m2 
  thick_m = thick_sum

  dens_m = (thick_m1*dens_m1 + thick_m2*dens_m2) / thick_sum
  ttime_sum = thick_m1/vp_m1 + thick_m2/vp_m2
  vp_m   = thick_sum / ttime_sum

  ttime_sum = thick_m1/vs_m1 + thick_m2/vs_m2
  vs_m   = thick_sum / ttime_sum
     
  return
  end
