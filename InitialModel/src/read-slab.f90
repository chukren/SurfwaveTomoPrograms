!=======================================================================
!  Read slab depth from cas_slab1.0_clip.xyz (USGS 2012)
!  
!  Nov 24 2014  chukren  Princeton
!=======================================================================
  module slabpara
  implicit none

  private 

  real,     parameter :: tiny = 1.0E-06
  integer,  parameter :: nlatmax = 651  ! fixed
  integer,  parameter :: nlonmax = 401  ! fixed
  real,     parameter :: slablonmin = -128.50
  real,     parameter :: slablonmax = -120.50
  real,     parameter :: slablatmin = 39.0
  real,     parameter :: slablatmax = 52.0
  real,     parameter :: slabxdim = 0.02
  real,     parameter :: slabydim = 0.02

  real,     dimension(1:nlatmax) :: slablat=0.
  real,     dimension(1:nlonmax) :: slablon=0.
  real,     dimension(1:nlonmax,1:nlatmax) :: hslab=-999.

  character(len = 70),parameter :: fslab= "cas_slab1.0_clip.xyz"

  public :: read_slab
  public :: read_slab_xyz

  contains

  !program readslabdep
  !  use slabpara
  !  implicit none
  !  real :: slon = -122.230
  !  real :: slat =   39.410
  !  real :: slabdep 
  !  call read_slab_xyz()
  !  call read_slab(slon,slat,slabdep)
  !  write(*,*)'slab depth =',slabdep
  !end

  subroutine read_slab(slon,slat,slabdep)
    !use slabpara
    implicit none

    integer :: i,j
    real, parameter :: slabnan = -999.
    real :: slabdep
    real :: slon, slat
    real :: londif, latdif

    if(slon > 180) slon = slon - 360

    londif = slon-slablonmin
    latdif = slat-slablatmin

    i = nint(londif/slabxdim) + 1
    j = nint(latdif/slabydim) + 1

    if(i > nlonmax .or. i < 0 .or. j > nlatmax .or. j < 0) then
      slabdep = slabnan
    else 
      slabdep = hslab(i,j)
    endif

    return
  end subroutine

  subroutine read_slab_xyz()
    !use slabpara
    !USE,INTRINSIC :: IEEE_ARITHMETIC ! for NaN check
    use, intrinsic :: ieee_arithmetic ! for NaN check
    implicit none
    real, parameter :: slabnan = -999.
    integer :: stat, UnitNum
    real :: xlon, xlat, xslab
    real :: difflat, difflon
    integer :: i, j

    UnitNum = 99

    open(unit=UnitNum, file=fslab, status='old', action='read' )

    do
      read (UnitNum, *, iostat=stat) xlon, xlat, xslab
      if (stat /= 0) exit

      if(xlon.ge.slablonmin .and. xlon.le.slablonmax .and. & 
         xlat.ge.slablatmin .and. xlat.le.slablatmax) then
        difflon = xlon - slablonmin
        difflat = xlat - slablatmin
        
        i = nint(difflon/slabxdim) + 1
        j = nint(difflat/slabydim) + 1

        slablon(i) = xlon
        slablat(j) = xlat
        
        ! check NaN
        if(.not. ieee_is_nan(xslab)) then
            hslab(i,j) = -1.0 * xslab
        else
            hslab(i,j) = slabnan 
        endif
        !write(*,1002)i,j,xlat,xlon,hslab(i,j)
1002 format(I5,1x,I5,1x,F12.3,1x,F12.3,1x,F12.3)

      endif
    end do

    close (UnitNum)

    return
   end subroutine 

  end module
