!=======================================================================
! Read seafloor age from NOAA's global dataset
! 
! Y.Y. Ruan Brown 08/29/2014
!
! 06/18/2015 chukren Princeton
!   Enlarged the area and changed the parameters
!======================================================================= 
module globalage
  implicit none

  private 

  !real, parameter :: lonmin = -132.
  !real, parameter :: lonmax = -120.
  !real, parameter :: latmin = 38.
  !real, parameter :: latmax = 52.

  real, parameter :: lonmin = -136.
  real, parameter :: lonmax = -118.
  real, parameter :: latmin = 36.
  real, parameter :: latmax = 54.

  real, parameter :: loninc = 0.1
  real, parameter :: latinc = 0.1

  integer, parameter :: nlon = 181 
  integer, parameter :: nlat = 181 

  real :: clon(nlon), clat(nlat)
  real :: cage(nlon,nlat)
  

  public :: readglobalcrustage
  public :: crustageglb

  contains

  !===
  !   read in whole model database in the designed area
  !===
  subroutine readglobalcrustage()
  use ieee_arithmetic
  implicit none
  character (len=70) :: fileage 
  real :: xlon, xlat, xage
  integer :: i, j
  integer :: n, maxline 
  
  !fileage = 'global.age.JdF.xy'
  fileage = 'global.age.JdF.large.xy'
  open(21, file = fileage)
  
  maxline = nlon*nlat

  do n = 1, maxline 
    read(21,*)xlon, xlat, xage
    i = nint( (xlon - lonmin)/loninc ) + 1
    j = nint( (xlat - latmin)/latinc ) + 1
    clon(i)   = xlon
    clat(j)   = xlat

    if(ieee_is_nan(xage)) then 
      cage(i,j) = -999.0
    else 
      cage(i,j) = xage
    endif

  enddo 

  close(21)

  return 
  end subroutine


  !===
  !   read the model pointwise
  !===
  subroutine crustageglb(lon, lat, age)
  implicit none
  
  real, parameter :: eps = 1.0E-06
  real :: lon, lat, age
  real :: x, y
  integer :: i1, i2, j1, j2 
  integer :: i, j

  x = (lon - lonmin)/loninc 
  y = (lat - latmin)/latinc

  i1 = floor(x) + 1
  j1 = floor(y) + 1
  i2 = i1 + 1
  j2 = j1 + 1

  if(i2.gt.nlon .or. j2.gt.nlat) then
    age = -999.0
  else
    !=== 
    !    Bilinear interpolation 
    !===
    age =       cage(i1,j1)*(clon(i2) - lon)*(clat(j2) - lat)
    age = age + cage(i2,j1)*(lon - clon(i1))*(clat(j2) - lat)
    age = age + cage(i1,j2)*(clon(i2) - lon)*(lat - clat(j1))
    age = age + cage(i2,j2)*(lon - clon(i1))*(lat - clat(j1))
    age = age/(loninc*latinc)

    if(age.lt.0.0) age = -999.0
  endif

  return
  end subroutine

end module
