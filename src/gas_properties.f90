!==========================================================
! This module contains functions for computing gas
! properties.
!==========================================================

module gas_properties
  implicit none
  private
  public :: mu, k

  !---------------------------------------------------------
  ! Function for computing molecular viscosity using
  ! Sutherland's law.
  !---------------------------------------------------------
  pure function mu(T)
    implicit none
    double precision, intent(in) :: T
    double precision             :: mu
    mu = 1.45d0*T**(3.0d0/2.0d0)/(T+110.0d0)*1.0e-6
  end function mu

  !---------------------------------------------------------
  ! Function for computing the thermal conductivity using
  ! temperature, Prandtl number and the specific heat at
  ! constant pressure.
  !---------------------------------------------------------
  pure function k(T,Pr,cp)
    implicit none
    double precision, intent(in) :: T,Pr,cp
    double precision             :: k
    k = cp/Pr*mu(T)
  end function k

end module gas_properties
