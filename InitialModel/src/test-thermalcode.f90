! Generate shear wavespeed profile from a thermal model
! Originally written by Donald Forsyth
! Modified by Y.Y. Ruan 08/16／2015
!=======================================================

  program testthermodel
  use thermodel
  implicit none
  
  real :: tmax, age, temp, vs, vp, Q
  real :: depth, thick
  real :: seisfreq, Tm, eta, J1, J2, Q_b, Vs_b

  integer :: i


  seisfreq = 0.02 ! 50 sec
      
  call setupzeroage(zeroaget,zerocoeff,pottemp,adbatgrad,srate)
  tmax = zeroaget(201)

  !write(*,*)'input age (Ma):'
  read(*,*)age
 
  do i = 1, 201
      depth = (i-1)  ! km  
      call tempest(zeroaget,zerocoeff,age,depth,tmax,temp)      
      !call mccarthytakei2(depth,temp,vs,Q)
      call yamauchi_takei_2016(depth,temp,seisfreq,vs,Q,Tm,eta,J1,J2,Q_b,Vs_b)
      write(*,*) depth,zeroaget(i),temp,vs,Q
  enddo    
 
  !depth = 60.
  !thick = 0.1
  !call thermantle(age,depth,thick,vs,vp)
  !write(*,*) depth,thick,vs,vp

  end
