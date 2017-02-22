!=====================================
! Read water depth
! 
! Input:  
!       real :: lon, lat
! Output: 
!       real :: water depth (m)
!
! Source:
!      etopo1_ice_g_i2
!=====================================
      module topopara
      implicit none 

      private 

      integer, parameter :: ncols = 1801
      integer, parameter :: nrows = 1801
      double precision, parameter :: xdim = 0.016666666667
      double precision, parameter :: ydim = 0.016666666667
      real, parameter :: lon0 = -150.0   ! lon = [-150, -120]
      real, parameter :: lat0 =   60.0   ! lat = [60, 30]
      integer, parameter :: errwat = 9999

      character(len = *), parameter :: ftopoxyz = 'JdFtopo.xyz'
      character(len = *), parameter :: ftopobin = 'JdFtopo.bin'

      integer :: watdep(ncols,nrows)
      real    :: grdlon(ncols)
      real    :: grdlat(nrows) 

      public :: gentopobin
      public :: readwatdep2
      public :: readwatdep
      public :: readtopobin

      contains

      subroutine gentopobin()
      !use topopara
      implicit none

      integer :: stat, i, j
      integer :: topo
      real :: xlon, xlat

      open(10,file=ftopoxyz,status='old')
      open(11,file=ftopobin,status='unknown',form='binary')

      !do j = 1, nrows
      !  do i = 1, ncols
      !    read(10,*,iostat=stat) grdlon(i),grdlat(j),watdep(i,j) 
      !    write(11) watdep(i,j) 
      !    if (stat /= 0) then 
      !      write(*,*)'The End!'
      !      exit
      !    endif 
      !  enddo
      !enddo       
      
      do 
        read (10, *, iostat=stat) xlon, xlat, topo

        i = nint( (xlon-lon0)/xdim ) + 1
        j = nint( (lat0-xlat)/ydim ) + 1
        watdep(i,j) = topo
        write(11) watdep(i,j)

        if (stat /= 0) then 
          !write(*,*)'The End of file!'
          exit
        endif
      enddo

      close(10)
      close(11)
      write(*,*)'Binary table file ',ftopobin,' is generated!' 

      write(*,1001) grdlon(ncols),grdlat(nrows),watdep(1801,1801)
1001 format(f8.3,1x,f8.3,1x,I5)

      return 
      end subroutine 


      !===
      !   read watdep(i,j) as 
      !===
      subroutine readwatdep2(lon,lat,watdepth)
      !use topopara
      implicit none

      real :: lon, lat
      integer :: i,j,watdepth

      i = nint((lon - lon0) / xdim) + 1
      if(i > ncols .or. i < 1) then
        print *,'Out of topo range lon = [-150W, -120W]! Set to',errwat
        watdepth = errwat
        return
      endif

      j = nint((lat0 - lat) / ydim) + 1
      if(j > nrows .or. j < 1) then
        print *,'Out of topo range lat = [60N, 30N]! Set to',errwat
        watdepth = errwat
        return
      endif

      watdepth = watdep(i,j)
      
      !write(*,1001) lon,lat,watdepth
1001 format(f8.3,1x,f8.3,1x,I5)

      return
      end subroutine

      !===
      !   read binary topo file at (lon,lat)
      !===
      subroutine readwatdep(lon,lat,watdepth)
      !use topopara
      implicit none

      real :: lon, lat
      integer :: i,j,watdepth
    
      open(10,file=ftopobin,status='old',form='binary')

      do j = 1, nrows
        do i = 1, ncols
          read(10) watdep(i,j)
        enddo
      enddo

      i = nint((lon - lon0) / xdim) + 1
      if(i > ncols .or. i < 1) then
        print *,'Out of range lon = [-150W, -120W]! set to',errwat
        watdepth = errwat
        return
      endif

      j = nint((lat0 - lat) / ydim) + 1
      if(j > nrows .or. j < 1) then
        print *,'Out of range lat = [60N, 30N]!, set to',errwat
        watdepth = errwat
        return
      endif

      watdepth = watdep(i,j)

      if(watdepth.le.0) then 
        watdepth = -watdepth
      else
        watdepth = 0
      endif

      !write(*,1001) lon,lat,watdepth
1001 format(f8.3,1x,f8.3,1x,I5)

      return
      end subroutine


      subroutine readtopobin()
      !use topopara
      implicit none
      integer :: i, j

      open(10,file=ftopobin,status='old',form='binary')

      do j = 1, nrows
        do i = 1, ncols
          read(10) watdep(i,j)
        enddo
      enddo

      return
      end subroutine

      end module
