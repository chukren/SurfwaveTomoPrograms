! Generate shear wavespeed and attenuation profile from a thermal model
! (McCarthy et al., 2011)
! Originally written by Donald Forsyth
! Modified by Y.Y. Ruan 08/16/2015
!
! Y.Y. Ruan 12/19/2016 
!   Add Grain Boundary Sliding peak (GBS) (Olugboji et al., 2013) 
!   Add Arrhenius attenuation
!
! Y.Y. Ruan 01/05/2017
!   Add near solidus melting effects on relaxation (Yamauchi & Takei et al., 2016) 
!
! Y.Y. Ruan 01/26/2017
!   Add peritotite solidus (Hirschmann, 2000)
!
!=======================================================

  module thermodel
  implicit none  
    
  private

  !===
  ! constants
  !===

  real, parameter :: R = 8.314462            ! gas constant
  integer, parameter :: maxvert = 400
  real, parameter :: pi = 3.14159265359
  real, parameter :: pottemp = 1350. ! potential temperature (c)
  real, parameter :: adbatgrad = 0.4 ! adiabatic gradient in deg C/km
  real, parameter :: srate = 30.0    ! spreading half rate, mm/yr

  real :: zeroaget(maxvert), zerocoeff(maxvert)

  !=== 
  ! subroutines 
  !===

  public :: setupzeroage
  public :: tempest
  public :: linear_solidus
  public :: predotite_solidus_Hirschmann2000
  public :: predotite_wetsolidus_Hirschmann2000
  public :: hydrostatic_pressure
  public :: unrelaxed_shear_modulus
  public :: unrelaxed_shear_modulus_issak_1992
  public :: unrelaxed_shear_modulus_PM2003
  public :: unrelaxed_shear_modulus_FJ2005
  public :: unrelaxed_shear_modulus_JF2010
  public :: shear_viscosity
  public :: mccarthy_etal_2011
  public :: yamauchi_takei_2016
  public :: thermantle
  public :: grain_boundary_sliding_atten
  public :: arrhenius_atten
 
  !=== 
  ! variables 
  !===

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
  !    call mccarthy_etal_2011(depth,temp,Vs,Q)
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

      ! extra parameters for yamauchi_takei_2016()
      real :: seisfreq,Tm,eta,J1,J2,Q_b,Vs_b 

      integer :: nl = 100 ! integral points
      integer :: i
     
      seisfreq = 0.02 ! 50 s 

      tmax = zeroaget(201)
      thick_inc = thick / nl
      thick_hf = thick_inc / 2.0
      t_sum = 0.0

      do i = 1, nl
        depth_i = depth + (i-1)*thick_inc + thick_hf
        call tempest(zeroaget,zerocoeff,age,depth_i,tmax,temp)
        call yamauchi_takei_2016(depth_i,temp,seisfreq,vs,Q,Tm,eta,J1,J2,Q_b,Vs_b)
        !call mccarthy_etal_2011(depth_i,temp,vs,Q,temp)
        t_sum = t_sum + thick_inc/vs 
      end do

      vs = thick / t_sum
      vp = vs * vpvs_ratio

      return 
    end subroutine


    subroutine hydrostatic_pressure(depth, pressure)
      implicit none

      real :: depth, pressure, pressure_gpa

      ! crudely calculate pressure at given depth - ignore density diff in crust
      pressure = depth * 3.2 * 9.8 * 1.0E+6
      pressure_gpa = pressure / 1.0E9 

      return
    end subroutine


    subroutine linear_solidus(depth, temperature)
      implicit none
      ! dry solidus
      real, parameter :: solidus0 = 1100.       ! dry solidus temperature at zero pressure (deg C)
      real, parameter :: solgrad = 3.5          ! gradient of dry solidus (degC/km) 
      real :: depth, temperature                ! km, C

      temperature = solidus0 + solgrad * depth

      return
    end subroutine


    subroutine predotite_solidus_Hirschmann2000(pressure, temperature)
      implicit none
      real, parameter :: coeff_a = -5.14   ! C/GPa^2
      real, parameter :: coeff_b = 132.9   ! C/GPa
      real, parameter :: coeff_c = 1120.7  ! C

      real :: pressure, temperature

      Temperature = coeff_c + coeff_b * pressure + coeff_a * pressure**2

      return
    end subroutine


    subroutine predotite_wetsolidus_Hirschmann2000(pressure, temperature)
      implicit none
      real, parameter :: coeff_a = -5.14   ! C/GPa^2
      real, parameter :: coeff_b = 101.9   ! C/GPa
      real, parameter :: coeff_c = 1120.7  ! C

      real :: pressure, temperature

      Temperature = coeff_c + coeff_b * pressure + coeff_a * pressure**2

      return
    end subroutine


    subroutine unrelaxed_shear_modulus_issak_1992(temperature, pressure, shear_modulus)
      implicit none

      !=== From Issak, JGR, 1992 and used in McCarthy et al. (2011)
      real, parameter :: Gur   = 82.45E+09      ! unrelaxed modulus (Pa)
      real, parameter :: dGudT = -13.6E+06      ! temp derivative Pa/K  [-1.38E-2, -1.36E-2]GPa/K
      real, parameter :: dGudP = 1.8            ! pressure derivative Pa/Pa

      real :: temperature, pressure
      real :: shear_modulus

      if (temperature .lt. 0 .or. temperature .gt. 6000) then
          write(*,*)"Error: Impossible temperature!"
          stop
      endif

      ! calculate unrelaxed modulus at given temperature and pressure
      shear_modulus = Gur + (temperature + 273.) * dGudT + pressure * dGudP

      return
    end subroutine


    subroutine unrelaxed_shear_modulus(temperature, pressure, shear_modulus)
      implicit none

      !=== Inverted from Priestley & McKenzie, EPSL, (2013) by Yamauchi & Takei
      ! JGR (2016)
      real, parameter :: Gur   = 72.45E+09      ! unrelaxed modulus (Pa)
      real, parameter :: dGudT = -10.94E+06     ! temp derivative Pa/degK
      !real, parameter :: dGudP = 1.987          ! pressure derivative Pa/Pa
      real, parameter :: dGudP = 1.5           ! pressure derivative Pa/Pa

      real :: temperature, pressure
      real :: shear_modulus     ! Pa

      if (temperature .lt. 0 .or. temperature .gt. 6000) then
          write(*,*)"Error: Impossible temperature!"
          stop
      endif

      ! calculate unrelaxed modulus at given temperature and pressure
      shear_modulus = Gur + (temperature + 273.) * dGudT + pressure * dGudP

      return
    end subroutine

    subroutine unrelaxed_shear_modulus_PM2003(temperature, pressure, shear_modulus)
      implicit none

      !=== Inverted from Priestley & McKenzie, EPSL, (2013)
      real, parameter :: Gur   = 72.66E+09      ! unrelaxed modulus (Pa)
      real, parameter :: dGudT = -8.71E+06     ! temp derivative Pa/degK
      real, parameter :: dGudP = 2.04          ! pressure derivative Pa/Pa

      real :: temperature, pressure
      real :: shear_modulus     ! Pa

      if (temperature .lt. 0 .or. temperature .gt. 6000) then
          write(*,*)"Error: Impossible temperature!"
          stop
      endif

      ! calculate unrelaxed modulus at given temperature and pressure
      shear_modulus = Gur + (temperature + 273.) * dGudT + pressure * dGudP

      return
    end subroutine


    subroutine unrelaxed_shear_modulus_JF2010(temperature, pressure, shear_modulus)
      implicit none

      !=== Inverted from Jackson & Fual, PEPI, (2010)
      real, parameter :: Gur   = 66.65E+09      ! unrelaxed modulus (Pa)
      real, parameter :: dGudT = -13.6E+06     ! temp derivative Pa/degK
      real, parameter :: dGudP = 1.8          ! pressure derivative Pa/Pa

      real :: temperature, pressure
      real :: shear_modulus     ! Pa

      if (temperature .lt. 0 .or. temperature .gt. 6000) then
          write(*,*)"Error: Impossible temperature!"
          stop
      endif

      ! calculate unrelaxed modulus at given temperature and pressure
      shear_modulus = Gur + (temperature + 273. - 1173) * dGudT + (pressure - 0.2E09) * dGudP

      return
    end subroutine


    subroutine unrelaxed_shear_modulus_FJ2005(temperature, pressure, grainsize, shear_modulus)
        implicit none
        ! From Faul & Jackson, 2005; Tabel 1
        ! Ju = Ju(P) + dln(Ju)

        real, parameter :: JuP = 0.0149             ! 1/GPa
        real, parameter :: dlnJudT = 9.1E-04        ! 1/K
        real, parameter :: dGudP = 1.8              
        real, parameter :: grainsize_ref = 1.0E-05  ! meter
        real, parameter :: grainsize_coef = -0.16   ! 
        real, parameter :: Tr = 1223                ! k

        real :: temperature, grainsize, Ju, pressure, shear_modulus
        real :: grainsize_term

        grainsize_term = (grainsize / grainsize_ref)**grainsize_coef

        Ju = JuP * (1.0 + dlnJudT * grainsize_term * (temperature - Tr))

        shear_modulus = 1.0E09 / Ju + dGudP * pressure

        return
    end subroutine

    subroutine shear_viscosity(temperature, pressure, premelt_coeff, viscosity)
      implicit none

      !=== Inverted from Priestley & McKenzie, EPSL, (2013) by Yamauchi & Takei
      ! JGR (2016)
      real, parameter :: Tr = 1200. + 273.      ! reference temperature
      real, parameter :: Pr = 1.5E+09           ! reference pressure  
      real, parameter :: V  = 7.913E-06         ! activation volume (m^3/mol)
      real, parameter :: E1 = 462.5             ! activation energy in kJ/mol  
      real, parameter :: eta0 = 6.22E+21        ! reference viscosity

      !=== From McCarthy et al. (2011) P.13
      !real, parameter :: Tr = 1200. + 273.     ! reference temperature
      !real, parameter :: Pr = 2.0E+09          ! reference pressure (???) 
      !real, parameter :: V  = 12.E-06          ! activation volume (m^3/mol)
      !real, parameter :: E1 = 505              ! activation energy in kJ/mol  
      !real, parameter :: eta0 = 6.6E+19        ! reference viscosity (6.6E19)

      !=== From Donald Forsyth (personal communication)
      !real, parameter :: Tr = 1425. + 273.     ! reference temperature (?)
      !real, parameter :: Pr = 2.0E+09          ! reference pressure  
      !real, parameter :: V  = 12.0E-06         ! activation volume (m^3/mol)
      !real, parameter :: E1 = 400              ! activation energy in kJ/mol  
      !real, parameter :: eta0 = 1.0E+21        ! reference viscosity (6.6E19)

      !=== lower down eta0 will reduce the vel and atten, linear relation

      real, parameter :: E  = E1 * 1000.        ! activation energy in J/mol
      real, parameter :: viscosity_max = 1.0E+23

      !=== Input/Output parameters
      real :: temperature, pressure, premelt_coeff, viscosity

      real :: energy_term, volume_term
      real :: temperature_c, eta_nomelt

      !=== temperature C to K 
      temperature_c = temperature + 273. 

      energy_term = exp((1. / temperature_c - 1./Tr) * E / R)
      volume_term = exp((pressure / temperature_c - Pr / Tr) * V / R)

      ! calculate viscosity (keeping pressure at reference pressure) eqn. (17)
      eta_nomelt = eta0 * energy_term * volume_term

      if(eta_nomelt .gt. viscosity_max) eta_nomelt = viscosity_max

      viscosity = eta_nomelt * premelt_coeff

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
        temp = zeroaget(idepth) + (zeroaget(idepth+1) - zeroaget(idepth)) * (depth + 1. - depthi)
      else
        tpikappal2 = -kappa * age * 3.1557E+07 * (pi/bot)**2
        sumt = 0.0
        do n = 1, 100
          sumt = sumt + zerocoeff(n) * sin(n* pi * depth / bot) * exp(tpikappal2 * n * n)
        enddo
        temp = sumt + tmax * depth/bot
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
    ! DWF (original) YYR (modified)
    !===========================================================
    subroutine mccarthy_etal_2011(z,T,Vs,Q,eta)
      implicit none

      real, parameter :: R= 8.314462            ! gas constant
      real, parameter :: twopi = 2.0*pi 

      real, parameter :: Gur   = 82.0E+09       ! unrelaxed modulus (Pa)
      real, parameter :: dGudT = -13.6E+06      ! temp derivative Pa/degK
      real, parameter :: dGudP = 1.8            ! pressure derivative Pa/Pa

      real, parameter :: Tr = 1425. + 273.      ! reference temperature
      real, parameter :: Pr = 2.0E+09           ! reference pressure  
      real, parameter :: V  = 12.0E-06          ! activation volume (m^3/mol)
      real, parameter :: E1 = 400.              ! activation energy in kJ/mol  
      real, parameter :: E  = E1*1000.          ! activation energy in J/mol
      real, parameter :: eta0 = 1.0E+21         ! reference viscosity (6.6E19)

      real, parameter :: solidus0 = 1100.       ! dry solidus temperature at zero pressure (deg C)
      real, parameter :: solgrad = 3.5          ! gradient of dry solidus (degC/km)  

      real :: a0,a1,a2,a3,a4,a5,a6
      real :: P, z, T
      real :: eta, Gu, dGup, taum, G1, taufac, xn, qinv, q, Vs
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
      ! Note the original code doesn't take into account the pressure
      ! effects by mistake. 
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

    !===========================================================      
    ! This program predicts shear velocity using results from 
    ! Yamauchi & Takei, JGR (2016), equation 11-18 and constants 
    ! in Table 4 and Fig 20.
    !===========================================================
    subroutine yamauchi_takei_2016(z,T,seisfreq,Vs,Q,Tm,eta,J1,J2,Q_b,Vs_b)
      implicit none

      !=== parameter are from Table 4 and Fig. 20
      real, parameter :: R= 8.314462            ! gas constant
      real, parameter :: twopi = 2.0*pi         !
      real, parameter :: Tr = 1200. + 273.      ! reference temperature
      real, parameter :: Pr = 1.5E+09           ! reference pressure  
      real, parameter :: V  = 7.913E-06         ! activation volume (m^3/mol)
      real, parameter :: E1 = 462.5             ! activation energy in kJ/mol  
      real, parameter :: E  = E1*1000.          ! activation energy in J/mol

      ! lower down the vel and atten, linear relation
      real, parameter :: eta0 = 6.22E+21        ! reference viscosity
    
      ! dry solidus
      real, parameter :: solidus0 = 1100.       ! dry solidus temperature at zero pressure (deg C)
      real, parameter :: solgrad = 3.5          ! gradient of dry solidus (degC/km) 

      real, parameter :: Tn_eta = 0.94          ! normalized temperature above which
      real, parameter :: gama = 5.0             ! the factor of extra reduction in eta
      real, parameter :: lambda = 26            ! melt fraction dependence 
      real, parameter :: density_m = 3.3E+03    ! mantle density used for shear wavespeed
      real, parameter :: tau_n_p = 6.0E-05      ! normalized time scale
      real, parameter :: amp_b = 0.664
      
      ! change the minimum of atten linearly
      real, parameter :: alpha = 0.38
      
      real :: a0,a1,a2,a3,a4,a5,a6
      real :: P, z, T, P_gpa
      real :: eta, Gu, dGup, tau_m, G1, taufac, xn, qinv, q, Vs
      real :: seisfreq, fn, fnln, tau_n
      real :: J1, J2, J1_b, J2_b
      real :: q_b, qinv_b, vs_b
      real :: Tn, Tm
      real :: A_eta ! exrta term in eqn. (17) to account for reduction of eta just below solidus 
      real :: amp_p, sigma_p ! amplitude and width of peak
      real :: melt_frac
      real :: T_melt, P_melt, eta_melt

      ! crudely calculate pressure at given depth - ignore density diff in crust
      P = z * 3.2 * 9.8 * 1.0E+6
      P_gpa = P / 1.0E9 ! GPa

      !call predotite_solidus_Hirschmann2000(P_gpa, Tm)
      call predotite_wetsolidus_Hirschmann2000(P_gpa, Tm)

      Tn = T / Tm

      !=== melt fraction
      !if (Tn < 1.0) then
      !  melt_frac = 0.0
      !else
      !  melt_frac = (Tn - 1.0) * .3
      !endif

      melt_frac = 0.0

      !=== debug 
      !write(*,*) z,"melt_frac:", melt_frac


      !=== amplitude of peak: amp_p ===
      if (Tn < 0.91) then
          amp_p = 0.01 
      else if (Tn < 0.96 .and. Tn >= 0.91) then
          amp_p = 0.01 + 0.4 * (Tn - 0.91)
      else if (Tn < 1.0 .and. Tn >= 0.96) then
          amp_p = 0.03
      else ! Tn > 1.0
          amp_p = 0.03 
      endif

      !=== width of peak: sigma_p ===
      if (Tn < 0.92) then
          sigma_p = 4.0 
      else if (Tn < 1.0 .and. Tn >= 0.92) then
          sigma_p = 4.0 + 37.5 * (Tn - 0.92)
      else
          sigma_p = 7.0
      end if
      
      !write(*,1009) "depth/Tn/melt/Ap/sigma:", z, Tn, melt_frac, amp_p, sigma_p
1009 format(a, f9.2, f6.2, f9.6, f7.3, f7.3)

      !=== extra term ===
      if (Tn < Tn_eta) then
          A_eta = 1.0
      else if (Tn < 1 .and. Tn >= Tn_eta) then
          A_eta = exp(-1.0 * (Tn - Tn_eta) / (Tn - Tn*Tn_eta) * log(gama))
      else
          A_eta = 1. / gama * exp(-lambda * melt_frac)
      end if
 
      ! calculate viscosity (keeping pressure at reference pressure) eqn. (17)
      call shear_viscosity(T, P, A_eta, eta)

      ! calculate unrelaxed shear modulus (mu_u) at temperature and pressure
      ! Temperature: K
      ! Pressure: GPa
      ! Gu: Gpa
      call unrelaxed_shear_modulus(T, P, Gu)
      !call unrelaxed_shear_modulus_PM2003(T, P, Gu)
      !call unrelaxed_shear_modulus_JF2010(T, P, Gu)
      !call unrelaxed_shear_modulus_issak_1992(T, P, Gu)
 
      ! grainsize set to 0.005 m
      !call unrelaxed_shear_modulus_FJ2005(T, P, 0.005, Gu)

      ! calculate characteristic Maxwell relaxation time
      tau_m = eta / Gu

      ! using 50 s as typical period, calculate normalized frequency
      fn = seisfreq * tau_m
      tau_n = 1. / (twopi * fn)

      fnln = log(tau_n_p / tau_n)


      J1 = 0.5 * sqrt(twopi)*amp_p*sigma_p*(1.0 - erf(fnln/sigma_p/sqrt(2.0)))
      J1 = (1.0 / Gu) * (1.0 + amp_b * tau_n**alpha / alpha + J1)
      J1_b = (1.0 / Gu) * (1.0 + amp_b * tau_n**alpha)

      J2 = 0.5*pi*(amp_b*tau_n**alpha + amp_p*exp(-1.0*fnln**2/(2.*sigma_p**2))) + tau_n
      J2 = (1.0 / Gu) * J2
      J2_b = (1.0 / Gu) * (0.5*pi*(amp_b*tau_n**alpha) + tau_n)

      if (fnln.gt.1.0E13) then
        J1 = 1.0 / Gu
      endif

      qinv = J2 / J1
      vs = sqrt(1.0/(J1 * density_m)) / 1000. ! convert to km/s

      qinv_b = J2_b / J1_b
      vs_b = sqrt(1.0/(J1_b * density_m)) / 1000.

      if (qinv .lt. 0.002) qinv = .002
      q = 1.0 / qinv

      if (qinv_b .lt. 0.002) qinv_b = 0.002
      q_b = 1.0 / qinv_b

      return
    end subroutine


    subroutine grain_boundary_sliding_atten(depth, temperature, dvv_gbs, Q)
      implicit none

      ! From Olugboji et al. (2013)
      real, parameter :: gas_constant = 8.314
      real, parameter :: activation_p = 10.0E-06  ! m^3/mol 
      real, parameter :: activation_e = 3.56E+05  ! 356 KJ/mol
      real, parameter :: grainsize_d0 = 13.4E-06  ! m
      real, parameter :: watercontent_cw0 = 1.0E-4 ! wt%
      real, parameter :: Tr = 1173. + 273.  ! reference temperature (K)
      real, parameter :: Pr = 0.2E+09       ! reference pressure (Pa) 
      !real, parameter :: tau_pr = 1 ! reference relaxiation time 
      real, parameter :: log_tau_pr = -3.21
      
      ! user defined parameter
      real :: grainsize_exponent = 1.4 ! 1.4 +/- 0.47
      real :: grainsize = 1.0E-3 ! 1 cm

      real :: relax_strength_dp = 0.0667
      real :: omega_seis = 1.0 / 50. ! seismic frequency 

      real :: tau_gbs, log_tau_gbs, omega_gbs, log_omega_gbs
      real :: gs_term, log_gs_term
      real :: wt_term, log_wt_term
      real :: watercontent
      real :: watercontent_exponent ! 0 dry, 1 wet (weak), 2 wet (strong)
      real :: ref_term, syn_term
      real :: freqratio, freqratio_sqr
      real :: P, depth, temperature
      real :: Qinv, Q
      real :: dvv, dvv_gbs

      ! crudely calculate pressure given depth - ignore density diff in crust
      P = depth*3.2*9.8*1.0E+6

      if (depth < 70.) then
          watercontent = 0.0001
          watercontent_exponent = 0
      else if (depth >= 70. .and. depth <= 115.) then
          watercontent = 0.002 + (0.0001 - 0.002) * (depth - 115.) / (50.0 - 115.0)
          watercontent_exponent = 1
      else 
          watercontent = 0.002
          watercontent_exponent = 0
      endif

      gs_term = (grainsize / grainsize_d0) ** (-1.d0 * grainsize_exponent)
      wt_term = (watercontent / watercontent_cw0) ** watercontent_exponent

      ! calculate relaxation time due to grain boundary sliding (gbs)
      ref_term = exp( (activation_e + Pr * activation_p) / (gas_constant * Tr) )
      syn_term = exp( (activation_e + P  * activation_p) / (gas_constant * (temperature + 273)))

      ! Olugboji (eqn 10)
      tau_gbs = 10 ** log_tau_pr + gs_term * wt_term * ref_term / syn_term
      
      omega_gbs = 1.0 / tau_gbs
      log_omega_gbs = log10(omega_gbs)

      freqratio = omega_seis * tau_gbs
      freqratio_sqr = freqratio ** 2

      Qinv = relax_strength_dp * freqratio / (1 + freqratio_sqr)
      Q  = 1. / Qinv

      dvv_gbs = 0.5d0 * relax_strength_dp / (1 + freqratio_sqr)

    end subroutine


    subroutine arrhenius_atten(depth, temperature, Q)
      implicit none
        
      real, parameter :: gas_constant = 8.314
      real, parameter :: a_factor = 9.0E-04      !  3.0E-04

      real, parameter :: activation_e = 4.0E+05  !  450 KJ/mol  
      real, parameter :: activation_v = 20.E-06  !  13.E-06 m^3/mol

      real :: omega, enthalpy, temperature, pressure
      real :: alpha, Q, depth

      alpha = 0.27
      omega = 2 * 3.14159265

      ! crudely calculate pressure given depth - ignore density diff in crust
      pressure = depth*3.2*9.8*1.0E+6

      enthalpy = activation_e + pressure * activation_v
 
      Q = a_factor * (omega**alpha) * exp(alpha * enthalpy / gas_constant / temperature)

      return
    end subroutine

  end module  
