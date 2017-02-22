! Generate shear wavespeed profile from a thermal model
! Originally written by Donald Forsyth
! Modified by Y.Y. Ruan 08/16/2015
!=======================================================

  module thermodel
  implicit none  
    
  private

  integer, parameter :: maxvert = 400
  real, parameter :: pi = 3.14159265359
  real, parameter :: pottemp = 1350. ! potential temperature (c)
  real, parameter :: adbatgrad = 0.4 ! adiabatic gradient in deg C/km
  real, parameter :: srate = 30.0    ! spreading half rate, mm/yr

  real :: zeroaget(maxvert), zerocoeff(maxvert)

  !=== subroutines ===
  public :: setupzeroage
  public :: tempest
  public :: mccarthytakei2
  public :: thermantle
 
  !=== variables ===
  public :: zeroaget, zerocoeff
  public :: pottemp, adbatgrad, srate

  !==== original code keep for info ===================================
      
  ! pottemp is potential temperature deg C, adbatgrad is adiabatic gradient in 
  ! deg C/km, srate is spreading half rate, mm/yr
  !pottemp = 1350.
  !adbatgrad = 0.4
  !srate = 30.0
      
  !call setupzeroage(zeroaget,zerocoeff,pottemp,adbatgrad,srate)
  !tmax = zeroaget(201)

  ! age in Ma, depth in km 
  !age = 0.0
  !do i = 1,200
  !    depth = (i-1)    
  !    call tempest(zeroaget,zerocoeff,age,depth,tmax,temp)      
  !    call mccarthytakei2(depth,temp,Vs,Q)
  !    write(*,*) depth,zeroaget(i),temp,Vs,Q
  !enddo    
 
  !end
  !===================================================================


  contains


    subroutine thermantle(age,depth,thick,vs,vp)
      implicit none
      
      real, parameter :: vpvs_ratio = 1.73
      real :: tmax, age, temp
      real :: depth, thick, vs, vp, Q
      real :: dep_inc, dep_hf, depth_i
      real :: t_sum, thick_inc, thick_hf
      integer :: nl = 100 ! integral points
      integer :: i
     
      tmax = zeroaget(201)
      thick_inc = thick / nl
      t_sum = 0.0

      do i = 1, nl
        depth_i = depth + thick_hf + (i-1)*thick_inc
        call tempest(zeroaget,zerocoeff,age,depth_i,tmax,temp)
        call mccarthytakei2(depth_i,temp,vs,Q)
        t_sum = t_sum + thick_inc/vs 
      end do

      vs = thick / t_sum
      vp = vs * vpvs_ratio

      return 
    end subroutine

    subroutine tempest(zeroaget,zerocoeff,age,depth,tmax,temp)
      !parameter (maxvert = 400)
      implicit none

      real, parameter :: kappa = 1.0E-06
      real, parameter :: bot = 200.
      !real, parameter :: pi = 3.14159
      
      real :: zeroaget(maxvert),zerocoeff(maxvert)
      real :: age, depth, tmax, temp, sumt 
      real :: tpikappal2  
      
      integer :: idepth, depthi
      integer :: n

      if (depth.gt.bot) write(*,*) depth, ' depth too great'

      !  if age lt 0.1 Ma, then get temp from zeroage profile, otherwise 
      !  use Fourier series

      if (age.lt.0.1) then
        idepth = depth+1.0
        depthi = idepth
        !  this interpolation assumes depths at 1 km intervals
        temp = zeroaget(idepth) + (zeroaget(idepth+1)-zeroaget(idepth))*(depth+1.-depthi)
      else
        tpikappal2 = -kappa*age*3.1557E+07*(pi/bot)**2
      
        sumt = 0.0
        do n = 1, 100
          sumt = sumt + zerocoeff(n)*sin(n*pi*depth/bot)*exp(tpikappal2*n*n)
        enddo
        temp = sumt + tmax*depth/bot
      endif
      return
    end subroutine
     
        
    subroutine setupzeroage(zeroaget,zerocoeff,pottemp,adbatgrad,srate)
      implicit none

      real :: zeroaget(maxvert), zerocoeff(maxvert)
      real :: tzeroage(maxvert)
      real :: pottemp,adbatgrad,srate
      real :: kappa,bot

      real :: upwell,d,z,sumprod
      real :: zavg, tavg
      integer :: i, n
      
      ! upwelling rate for flat plate separation beneath axis is 2*srate/pi
      upwell = -(2.0*srate/pi)/3.15576E+10

      ! kappa is thermal diffusivity
      kappa = 1.0E-6

      d = upwell/kappa

      ! calculate temperature every 1 km from 0 to 200 km
      ! but make top 3 km zero to crudely mimic hydrothermal circulation cooling
      zeroaget(1) = 0.0
      zeroaget(2) = 0.0
      zeroaget(3) = 0.0       

      do i = 4, 201
        z = (i-4)*1000.
        zeroaget(i) = pottemp*(1.0 - exp(d*z))+adbatgrad*(i-1)
        ! write(*,*) z, zeroaget(i)
      enddo

      ! subtract steady-state temperature to get transient temp at axis, assuming
      ! temp fixed at bottom
      bot = 200.

      do i = 1, 201
        z = (i-1)
        tzeroage(i) = zeroaget(i) - zeroaget(201)*z/bot
        ! write(*,*) z, zeroaget(i),tzeroage(i)
      enddo

      ! find first 50 coefficients of Fourier sine expansion of transient temp
      ! would be faster to do FFT
      do n = 1, 100
         sumprod = 0.0
         do i = 1, 200 
           zavg = i - 0.5
           tavg = 0.5*(tzeroage(i)+tzeroage(i+1))
           sumprod = sumprod+ tavg*sin(n*pi*zavg/bot)
         enddo
         zerocoeff(n) = sumprod*2.0/bot
         !write(*,*) n, zerocoeff(n)
      enddo

      return
    end subroutine


    !===========================================================      
    ! This program (McCarthyTakeiAttn.f) predicts shear velocity 
    ! using universal curve from McCarthy et al, JGR 2011  equation 
    ! 26 and constants for olivine.
    ! parameters are temperature, viscosity at reference temperature, 
    ! ref temp, and activation energy
    !===========================================================
    subroutine mccarthytakei2(z,T,Vs,Q)
      implicit none

      real, parameter :: R= 8.314462            ! gas constant
      real, parameter :: twopi = 2.0*pi 
      real, parameter :: Gur   = 82.0E+09       ! unrelaxed modulus
      real, parameter :: dGudT = -13.6E+06      ! temp derivative Pa/degK
      real, parameter :: dGudP = 1.8            ! pressure derivative Pa/Pa
      real, parameter :: Tr = 1425. + 273.      ! reference temperature
      real, parameter :: Pr = 2.0E+09           ! reference pressure  
      real, parameter :: V  = 12.0E-06          ! activation volume
      real, parameter :: E1 = 400.              ! activation energy in kJ/mol  
      real, parameter :: E  = E1*1000. 
      real, parameter :: eta0 = 1.0E+21         ! reference viscosity
      real, parameter :: solidus0 = 1100.       ! dry solidus temperature at zero pressure (deg C)
      real, parameter :: solgrad = 3.5          ! gradient of dry solidus (degC/km)  

      real :: a0,a1,a2,a3,a4,a5,a6
      real :: P, z, T
      real :: eta, Gu, taum, G1, taufac, xn, qinv, q, Vs
      real :: fn, fnln, taun
      real :: j2

      ! coefficients for calculating modulus
      a0 =  0.55097
      a1 =  0.054332
      a2 = -0.0023615
      a3 = -5.7175E-05
      a4 =  9.9473E-06
      a5 = -3.4761E-07
      a6 =  3.9461E-09

      ! crudely calculate pressure given depth - ignore density diff in crust
      P = z*3.2*9.8*1.0E+6

      ! calculate viscosity (keeping pressure at reference pressure)
      eta = eta0*exp((E+P*V)/(R*(T+273)))/exp((E+Pr*V)/(R*Tr))

      ! reduce viscosity if melt present, i.e., depth > dry solidus
      if (T.gt.(solidus0+solgrad*z)) then
        eta = eta/100.
      endif

      ! calculate unrelaxed modulus at temperature and pressure
      ! Note in previous code the pressure term dGudP is not considered 
      Gu = Gur + (T+273.)*dGudT + P*dGudP
      ! write(*,*) Gu, Gur, T,P,dGudT,dGuP


      ! calculate characteristic Maxwell relaxation time
      taum = eta/Gu

      ! using 50 s as typical period, calculate normalized frequency
      fn = 0.02*taum
      taun = 1./(twopi*fn)
      fnln = alog(fn)

      ! use McCarthy's polynomial form to find relaxed modulus G1
      if (fnln.gt.1.0E13) then
        G1 = Gu
      else
        G1 = Gu*(a0 + fnln*a1 + a2*fnln**2 +a3*fnln**3 + a4*fnln**4 + a5*fnln**5 +a6*fnln**6)
      endif

      ! calculate Xn and Q using McCarthy's eqn. 25 & 18
      taufac = 0.39 - 0.28/(1.0+2.6*taun**0.1)
      xn = 0.32*taun**taufac
      j2 = (1./Gu)*(pi*xn/2. + 1./(2.*pi*fn))
      qinv = j2*Gu

      if (qinv.lt..002) qinv = .002
      q = 1./qinv
      Vs = sqrt(g1/3.3E03)/1000.
      !write(*,*) j2,Gu,q, Vs
      return
    end subroutine

  end module  
