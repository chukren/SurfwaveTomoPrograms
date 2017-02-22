!=======================================================================
!  Estimate sediment thickness for a station from sediment thickness 
!  according to NOAA.
!
!  Nov 19 2013  chukren  Brown
!=======================================================================
  module sedpara
  implicit none
  real,     parameter :: tiny = 1.0E-06
  integer,  parameter :: nmax = 361
  real,     parameter :: lonmin = -150.0
  real,     parameter :: lonmax = -120.0
  real,     parameter :: latmin = 30.0
  real,     parameter :: latmax = 60.0
  real,     parameter :: xdim = 8.35E-02
  real,     parameter :: ydim = 8.33E-02

  real,     dimension(1:nmax) :: lat=0.
  real,     dimension(1:nmax) :: lon=0.
  real,     dimension(1:nmax,1:nmax) :: hsed=9999.

  !character(len = 70),parameter :: fsed = "sedthick_values.xyz"
  character(len = 70),parameter :: fsed = "JdFsed.xyz"
  end module


  subroutine read_sed(slon,slat,sedthck)
    use sedpara
    implicit none

    integer, parameter :: nx = 4  ! nx * xdim (degree) 
    integer, parameter :: ny = 4  ! ny * ydim (degree)

    integer :: i,j,ix,iy,icount

    integer :: sedthck
    real :: slon, slat
    real :: londif, latdif, sedsum, diff
    real :: thick, thickmin, thickmax

    if(slon > 180) slon = slon - 360

    !call read_sed(10, fsed)

    thickmin = 20000.
    thickmax = 0.

    londif = slon-lonmin
    latdif = slat-latmin

    i = nint(londif/xdim) + 1
    j = nint(latdif/ydim) + 1


    diff = abs(hsed(i,j) - 9999.0)
    if(diff < tiny .or. i > nmax .or. j > nmax) then
      sedthck = 0
    else 
      !=== averaging ===
      icount = 0
      sedsum = 0.

      do ix = i-nx, i+nx
        if(ix < 1 .or. ix > nmax) cycle
        do iy = j-ny, j+ny
          if(iy < 1 .or. iy > nmax) cycle
          diff = abs(hsed(ix,iy) - 9999.0) 
          if(diff < tiny) cycle
          sedsum = sedsum + hsed(ix,iy)
          icount = icount + 1
        enddo
      enddo
      sedthck = nint(sedsum/icount)
    endif

    !write(*,1001)'sediment thickness =',sedthck
1001 format(a21,1x,I5)

    return
  end

  subroutine read_sed_xyz()
    use sedpara
    implicit none
    integer :: stat, UnitNum
    integer :: Numx, Numy
    real :: xlon, xlat, xsed
    real :: difflat, difflon
    integer :: i, j

    UnitNum = 10

    open(unit=UnitNum, file=fsed, status='old', action='read' )

    Numx = 0
    Numy = 0

    do
      read (UnitNum, *, iostat=stat) xlon, xlat, xsed 
      if (stat /= 0) exit

      if(xlon.ge.lonmin .and. xlon.le.lonmax .and. xlat.ge.latmin .and. xlat.le.latmax) then
        difflon = xlon - lonmin
        difflat = xlat - latmin

        i = nint(difflon/xdim) + 1
        j = nint(difflat/ydim) + 1

        Lon(i) = xlon
        Lat(j) = xlat
        hsed(i,j) = xsed

        !write(*,1002), i, j, xlon, xlat, xsed
1002 format(I5,1x,I5,1x,F12.5,1x,F12.5,1x,F8.2)

      endif
    end do

    close (UnitNum)

    return
   end subroutine 
