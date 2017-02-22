!========================================================================
!  Read the mantle model from 
!   ncdump -f f wUS-SH-2010_percent.nc > wUS-SH-2010.xyz
!   edit wUS-SH-2010.xyz for fortran readible (remove head , : // and ())
!  
!  This module contains one function
!    read_wUS_xyz()  
!  and model paremeters
!    lon, lat, dep, vp, vs 
!
!  Modified 06/03/2015 chukren Princeton
!========================================================================

  module wUSpara
  implicit none

  private

  !=== dimension of data ===
  integer, parameter :: ndim = 3
  integer, parameter :: ndep = 19 
  integer, parameter :: nlat = 93 
  integer, parameter :: nlon = 122
  
  character (len = *), parameter :: dep_name = "depth"
  character (len = *), parameter :: lat_name = "latitude"
  character (len = *), parameter :: lon_name = "longitude"

  real, parameter :: wus_lat_min = 27.50
  real, parameter :: wus_lat_max = 50.50
  real, parameter :: wus_lat_inc = 0.25
  real, parameter :: wus_lon_min = -125.75
  real, parameter :: wus_lon_max = -95.50
  real, parameter :: wus_lon_inc = 0.25
  real, parameter :: wus_dep_min = 45
  real, parameter :: wus_dep_max = 885

  !=== coordinates ===
  real :: wuslat(nlat)
  real :: wuslon(nlon) 
  real :: wusdep(ndep)

  !=== dvp and dvs ===
  character (len = *), parameter :: dvp_name = "dvp"
  character (len = *), parameter :: dvs_name = "dvs"

  real :: wusdvp(ndep,nlat,nlon)
  real :: wusdvs(ndep,nlat,nlon)
  
  character (len=*), parameter :: fname = "wUS-SH-2010.xyz" 

  !===
  !   Public the subroutine and model parameters
  !===
  public :: read_wUS_xyz
  public :: wusdvp, wusdvs
  public :: wuslat, wuslon, wusdep

  contains

  !===
  !   Read wUS model from a xyz file
  !=== 
  subroutine read_wUS_xyz()
  implicit none
  integer :: iidep, iilat, iilon
  integer :: i, j, k, stat
  integer :: nndep, nnlat, nnlon, ntotal
  integer :: unitnum = 30
  character :: dummy
 
  open(unit=unitnum,file=fname,status='old',action='read')

  do 
    !=== get the coordinates === 
    read(unitnum,*) nndep
    if(nndep .ne. ndep) then
      write(*,*)'depth layer wrong!'
      exit
    endif
    do iidep = 1, ndep
      read(unitnum,*) wusdep(iidep),dummy,i
    enddo
    
    read(unitnum,*) nnlat
    if(nnlat .ne. nlat) then
      write(*,*)'latitude dimension wrong!'
      exit
    endif
    do iilat = 1, nlat
      read(unitnum,*) wuslat(iilat),dummy,j
    enddo        
    
    read(unitnum,*) nnlon
    if(nnlon .ne. nlon) then
      write(*,*)'longitude dimension wrong!'
      exit
    endif
    do iilon = 1, nlon
      read(unitnum,*) wuslon(iilon),dummy,k
    enddo        
   
    !write(*,*)i,j,k  

    !=== get the data ===
    ! dvp(dep,lat,lon), dvs(dep,lat,lon)
    !====================
    ntotal = ndep * nlat * nlon    
    read(unitnum,*) ntotal
    do iidep = 1, ndep
      do iilat = 1, nlat
        do iilon = 1, nlon
          read(unitnum,*) wusdvp(iilon,iilat,iidep),dummy, k, j, i
        enddo
      enddo
    enddo

    !write(*,*)'dvp',k,j,i
    
    read(unitnum,*) ntotal
    do iidep = 1, ndep
      do iilat = 1, nlat
        do iilon = 1, nlon
          read(unitnum,*) wusdvs(iilon,iilat,iidep),dummy, k, j, i
        enddo
      enddo
    enddo
 
    !write(*,*)wusdvs(10,20,12),'10,20,12'

    read(unitnum,*,iostat=stat)
    if(stat /= 0) exit
  enddo
  
  close(unitnum)
  print *,"reading ",fname," done!"

  return
  end subroutine

  
  end module 
