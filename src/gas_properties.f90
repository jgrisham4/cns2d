!==========================================================
! This module contains functions for computing gas
! properties.
!
! Author: James Grisham
! Date: 01/13/2017
!==========================================================

module gas_properties
  implicit none
  private
  public :: mu, k

  ! Assuming constant Prandtl number
  !double precision, parameter :: Pr = 0.72d0

  ! Constants used in MMS
  double precision, parameter :: Pr = 1.0d0  ! MMS
  double precision, parameter :: mu_mms = 20.0d0

  contains

    !---------------------------------------------------------
    ! Function for computing molecular viscosity using
    ! Sutherland's law.
    !---------------------------------------------------------
    pure function mu(T)
      implicit none
      double precision, intent(in) :: T
      double precision             :: mu
      !mu = 1.45d0*T**(3.0d0/2.0d0)/(T+110.0d0)*1.0e-6
      mu = mu_mms
    end function mu

    !---------------------------------------------------------
    ! Function for computing the thermal conductivity using
    ! temperature, Prandtl number and the specific heat at
    ! constant pressure.
    !---------------------------------------------------------
    pure function k(T,g,R)
      implicit none
      double precision, intent(in) :: T,g,R
      double precision             :: k,cp
      cp = g*R/(g-1.0d0)
      k = mu(T)*cp/Pr
    end function k

end module gas_properties
