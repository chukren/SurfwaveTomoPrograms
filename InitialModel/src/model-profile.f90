! Test code to read a profile of 3D initial model
! chukren 08/18/2015
!
  program modelprofile
  implicit none 

  real, parameter :: eps = 1.0E-06
  real :: b_lat, e_lat, lat_inc
  real :: b_lon, e_lon, lon_inc
  real :: latx, lonx 
  real :: lat0, lon0
  real :: diff
  real :: thick, depth
  real :: rho, vp, vs 

  integer :: nlat, nlon, npts, nlay
  integer :: i, n 

  character(len=80) :: finp,fout

  write(*,*)'Input 3D model file:'
  read(*,*) finp

  write(*,*)'Input lat of profile'
  read(*,*) lat0

  write(*,*)'Output file for ploting:'
  read(*,*) fout

  open(11, file=finp, status='old', action='read')
  open(12, file=fout)

  read(11,*) b_lat, e_lat, lat_inc, b_lon, e_lon, lon_inc 

  nlat = nint((e_lat - b_lat)/lat_inc) + 1
  nlon = nint((e_lon - b_lon)/lon_inc) + 1
  npts = nlat * nlon

  do i = 1, npts

    read(11,*)lonx, latx    
    read(11,*)nlay

    diff = abs(latx - lat0)
    
    if(diff < eps) then
      ! generate depth point for plotting
      depth = 0.0
      do n = 1, nlay
        read(11,*) thick, rho, vp, vs
        depth = depth + thick
        write(12,1001) lonx, depth
1001 format(F10.3,1x,F7.3)     
      enddo
    else
      ! skip the point  
      do n = 1, nlay
        read(11,*)     
      enddo      
    endif 

  enddo

  close(11)
  close(12)
  end program 
