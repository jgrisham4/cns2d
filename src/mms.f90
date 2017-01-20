!==========================================================
! This module contains the source terms that are used
! in the method of manufactured solutions and the exact
! solutions.  The source terms are computed in Mathematica
! and written out to files in ../misc/source_*.  The exact
! solutions are also written out to files  in ../misc/
! from Mathematica.
!
! Author: James Grisham
! Date: 01-01-2017
!==========================================================

module mms
  implicit none
  private
  public :: rho_e,u_e,v_e,et_e,s_continuity,s_xmom,s_ymom,s_energy,dudx_e,dudy_e,dvdx_e,dvdy_e,dTdx_e,dTdy_e

  ! Parameters for MMS
  double precision, parameter :: rho0 = 0.5d0
  double precision, parameter :: u0   = 1.0d0
  double precision, parameter :: v0   = 0.1d0
  double precision, parameter :: et0  = 0.5d0
  double precision, parameter :: g    = 1.4d0
  double precision, parameter :: Rgas = 1.0d0
  double precision, parameter :: k    = 1.0d0
  double precision, parameter :: mu   = 0.3d0
  double precision, parameter :: eps  = 0.5d0

  contains

    !-------------------------------------------------------
    ! Function for rho(x,y)
    !-------------------------------------------------------
    pure function rho_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = rho0*(1.5 + Sin(x**2 + y**2))
    end function rho_e

    !-------------------------------------------------------
    ! Function for u(x,y)
    !-------------------------------------------------------
    pure function u_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = u0*(eps + Sin(x**2 + y**2))
    end function u_e

    !-------------------------------------------------------
    ! Function for v(x,y)
    !-------------------------------------------------------
    pure function v_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = v0*(eps + Cos(x**2 + y**2))
    end function v_e

    !-------------------------------------------------------
    ! Function for et(x,y)
    !-------------------------------------------------------
    pure function et_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = et0*(1.5 + Cos(x**2 + y**2))
    end function et_e

    !-------------------------------------------------------
    ! Function for computing du/dx
    !-------------------------------------------------------
    pure function dudx_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = 2*u0*x*Cos(x**2 + y**2)
    end function dudx_e

    !-------------------------------------------------------
    ! Function for computing du/dy
    !-------------------------------------------------------
    pure function dudy_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = 2*u0*y*Cos(x**2 + y**2)
    end function dudy_e

    !-------------------------------------------------------
    ! Function for computing dv/dx
    !-------------------------------------------------------
    pure function dvdx_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = -2*v0*x*Sin(x**2 + y**2)
    end function dvdx_e

    !-------------------------------------------------------
    ! Function for computing dv/dy
    !-------------------------------------------------------
    pure function dvdy_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = -2*v0*y*Sin(x**2 + y**2)
    end function dvdy_e

    !-------------------------------------------------------
    ! Function for computing dT/dx
    !-------------------------------------------------------
    pure function dTdx_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = ((-1 + g)*(-2*et0*x*Sin(x**2 + y**2) + &
        (4*v0**2*x*(eps + Cos(x**2 + y**2))*Sin(x**2 + y**2) - &
        4*u0**2*x*Cos(x**2 + y**2)*(eps + Sin(x**2 + y**2)))/2.))/Rgas
    end function dTdx_e

    !-------------------------------------------------------
    ! Function for computing dT/dy
    !-------------------------------------------------------
    pure function dTdy_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = ((-1 + g)*(-2*et0*y*Sin(x**2 + y**2) + &
        (4*v0**2*y*(eps + Cos(x**2 + y**2))*Sin(x**2 + y**2) - &
        4*u0**2*y*Cos(x**2 + y**2)*(eps + Sin(x**2 + y**2)))/2.))/Rgas
    end function dTdy_e

    !-------------------------------------------------------
    ! Function for computing the source term for the
    ! continuity equation
    !-------------------------------------------------------
    pure function s_continuity(x,y) result(s)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: s
      s = rho0*(v0*y*(2*Cos(2*(x**2 + y**2)) - 3*Sin(x**2 + y**2)) + &
        Cos(x**2 + y**2)*((3 + 2*eps)*u0*x + 2*eps*v0*y + &
        4*u0*x*Sin(x**2 + y**2)))
    end function s_continuity

    !-------------------------------------------------------
    ! Function for computing the source term for the
    ! x-momentum equation
    !-------------------------------------------------------
    pure function s_xmom(x,y) result(s)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: s

      ! This is the version for Lmv1
      s = (-3*(-1 + g)*rho0*v0**2*x*Cos(x**2 + y**2)**3 + &
        Sin(x**2 + y**2)*(-9*(-1 + g)*rho0*(et0 - eps*v0**2)*x - &
        4*mu*(v0 - 4*u0*x**2 + 2*u0*x*y - 3*u0*y**2) - &
        6*(-1 + g)*rho0*(et0 - eps*v0**2)*x*Sin(x**2 + y**2)) + &
        3*rho0*Cos(x**2 + y**2)**2*&
        (2*(-1 + g)*(et0 - eps*v0**2)*x + 3*u0*v0*y + &
        2*u0*v0*y*Sin(x**2 + y**2)) + &
        Cos(x**2 + y**2)*(3*eps*(6 + eps - (3 + eps)*g)*rho0*u0**2*x + &
        3*(-1 + g)*rho0*(3*et0 - eps**2*v0**2)*x + 4*mu*v0*(x - 2*y)*y + &
        u0*(-14*mu + 9*eps*rho0*v0*y) + &
        3*rho0*Sin(x**2 + y**2)*&
        ((6 + 6*eps - 3*g - 4*eps*g)*u0**2*x + 3*(-1 + g)*v0**2*x + &
        2*eps*u0*v0*y + ((5 - 3*g)*u0**2 + 2*(-1 + g)*v0**2)*x*&
        Sin(x**2 + y**2))))/3.

      ! This is the version for Lmv2
      !s = (-3*(-1 + g)*rho0*v0**2*x*Cos(x**2 + y**2)**3 + &
      !  3*rho0*Cos(x**2 + y**2)**2*&
      !  (2*et0*(-1 + g)*x + v0*(-2*eps*(-1 + g)*v0*x + (3 + 2*eps)*u0*y) + &
      !  4*u0*v0*y*Sin(x**2 + y**2)) + &
      !  Cos(x**2 + y**2)*(-3*eps*(3 + eps)*(-3 + g)*rho0*u0**2*x + &
      !  3*(-1 + g)*rho0*(3*et0 - eps**2*v0**2)*x + 4*mu*v0*(x - 2*y)*y + &
      !  u0*(-14*mu + 3*eps*(3 + 2*eps)*rho0*v0*y) + &
      !  3*rho0*Sin(x**2 + y**2)*&
      !  (-((3 + 4*eps)*(-3 + g)*u0**2*x) + 3*(-1 + g)*v0**2*x + &
      !  4*eps*u0*v0*y + (-3*(-3 + g)*u0**2 + 2*(-1 + g)*v0**2)*x*&
      !  Sin(x**2 + y**2))) + &
      !  Sin(x**2 + y**2)*(-9*rho0*(et0*(-1 + g)*x + &
      !  eps*v0*(v0*x - g*v0*x + u0*y)) - &
      !  4*mu*(v0 + u0*(-4*x**2 + 2*x*y - 3*y**2)) - &
      !  3*rho0*Sin(x**2 + y**2)*&
      !  (2*et0*(-1 + g)*x + v0*(-2*eps*(-1 + g)*v0*x + (3 + 2*eps)*u0*y) + &
      !  2*u0*v0*y*Sin(x**2 + y**2))))/3.
    end function s_xmom

    !-------------------------------------------------------
    ! Function for computing the source term for the
    ! y-momentum equation
    !-------------------------------------------------------
    pure function s_ymom(x,y) result(s)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: s

      ! This is the version for Lmv1
      s = (6*(-1 + g)*rho0*(et0 - eps*v0**2)*y*Cos(x**2 + y**2)**2 - &
        3*(-1 + g)*rho0*v0**2*y*Cos(x**2 + y**2)**3 + &
        Cos(x**2 + y**2)*(3*(-1 + g)*rho0*&
        (3*et0 - eps*((3 + eps)*u0**2 + eps*v0**2))*y + &
        4*mu*(u0 + v0*(3*x**2 - 2*x*y + 4*y**2)) + &
        3*rho0*y*Sin(x**2 + y**2)*&
        (-((3 + 4*eps)*(-1 + g)*u0**2) + 3*(-2 + g)*v0**2 + &
        (-3*(-1 + g)*u0**2 + 2*(-2 + g)*v0**2)*Sin(x**2 + y**2))) + &
        Sin(x**2 + y**2)*(2*mu*(7*v0 + 2*u0*x*(-2*x + y)) - &
        9*rho0*(et0*(-1 + g)*y + eps*v0*(u0*x - (-2 + g)*v0*y)) - &
        3*rho0*Sin(x**2 + y**2)*&
        ((3 + 2*eps)*u0*v0*x + 2*et0*(-1 + g)*y - 2*eps*(-2 + g)*v0**2*y + &
        2*u0*v0*x*Sin(x**2 + y**2))))/3.

      ! This is the version for Lmv2
      !s = (-3*(-3 + g)*rho0*v0**2*y*Cos(x**2 + y**2)**3 + &
      !  3*rho0*Cos(x**2 + y**2)**2*&
      !  ((3 + 2*eps)*u0*v0*x + 2*et0*(-1 + g)*y - 2*eps*(-3 + g)*v0**2*y + &
      !  4*u0*v0*x*Sin(x**2 + y**2)) + &
      !  Sin(x**2 + y**2)*(2*mu*(7*v0 + 2*u0*x*(-2*x + y)) - &
      !  9*rho0*(et0*(-1 + g)*y + eps*v0*(u0*x - (-3 + g)*v0*y)) - &
      !  3*rho0*Sin(x**2 + y**2)*&
      !  ((3 + 2*eps)*u0*v0*x + 2*et0*(-1 + g)*y - 2*eps*(-3 + g)*v0**2*y + &
      !  2*u0*v0*x*Sin(x**2 + y**2))) + &
      !  Cos(x**2 + y**2)*(3*rho0*(3*et0*(-1 + g)*y + &
      !  eps*((3 + 2*eps)*u0*v0*x - (3 + eps)*(-1 + g)*u0**2*y - &
      !  eps*(-3 + g)*v0**2*y)) + &
      !  4*mu*(u0 + v0*(3*x**2 - 2*x*y + 4*y**2)) + &
      !  3*rho0*Sin(x**2 + y**2)*&
      !  (4*eps*u0*v0*x - (3 + 4*eps)*(-1 + g)*u0**2*y + 3*(-3 + g)*v0**2*y + &
      !  (-3*(-1 + g)*u0**2 + 2*(-3 + g)*v0**2)*y*Sin(x**2 + y**2))))/3.
    end function s_ymom

    !-------------------------------------------------------
    ! Function for computing the source term for the
    ! energy equation
    !-------------------------------------------------------
    pure function s_energy(x,y) result(s)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: s
      s = ((3*et0*(Rgas*rho0*(9*u0*x + 6*eps*u0*x + v0*y + 6*eps*v0*y) + &
        8*(-1 + g)*k*(x**2 + y**2)) + &
        4*eps*(6*(-1 + g)*k*(u0**2 - v0**2*(x**2 + y**2)) + &
        mu*Rgas*(-7*u0**2 + 2*u0*v0*(1 + x*y - 2*y**2) + &
        2*v0**2*(3*x**2 - 2*x*y + 4*y**2))))*Cos(x**2 + y**2) + &
        2*(12*(-1 + g)*k*(u0**2 - v0**2)*(x**2 + y**2) + &
        Rgas*(u0*(4*mu*v0 + 3*(3 + 2*eps)*et0*rho0*x) - &
        4*mu*u0**2*(4*x**2 - 2*x*y + 3*y**2) + &
        v0*(3*(3 + 2*eps)*et0*rho0*y + 4*mu*v0*(3*x**2 - 2*x*y + 4*y**2))&
        ))*Cos(2*(x**2 + y**2)) + 9*et0*Rgas*rho0*v0*y*Cos(3*(x**2 + y**2)) - &
        24*et0*k*Sin(x**2 + y**2) + 24*et0*g*k*Sin(x**2 + y**2) - &
        8*eps*mu*Rgas*u0*v0*Sin(x**2 + y**2) + 24*eps*k*v0**2*Sin(x**2 + y**2) - &
        24*eps*g*k*v0**2*Sin(x**2 + y**2) + &
        28*eps*mu*Rgas*v0**2*Sin(x**2 + y**2) - &
        3*et0*Rgas*rho0*u0*x*Sin(x**2 + y**2) - &
        18*eps*et0*Rgas*rho0*u0*x*Sin(x**2 + y**2) + &
        24*eps*k*u0**2*x**2*Sin(x**2 + y**2) - &
        24*eps*g*k*u0**2*x**2*Sin(x**2 + y**2) + &
        32*eps*mu*Rgas*u0**2*x**2*Sin(x**2 + y**2) - &
        16*eps*mu*Rgas*u0*v0*x**2*Sin(x**2 + y**2) - &
        27*et0*Rgas*rho0*v0*y*Sin(x**2 + y**2) - &
        18*eps*et0*Rgas*rho0*v0*y*Sin(x**2 + y**2) - &
        16*eps*mu*Rgas*u0**2*x*y*Sin(x**2 + y**2) + &
        8*eps*mu*Rgas*u0*v0*x*y*Sin(x**2 + y**2) + &
        24*eps*k*u0**2*y**2*Sin(x**2 + y**2) - &
        24*eps*g*k*u0**2*y**2*Sin(x**2 + y**2) + &
        24*eps*mu*Rgas*u0**2*y**2*Sin(x**2 + y**2) - &
        12*k*u0**2*Sin(2*(x**2 + y**2)) + 12*g*k*u0**2*Sin(2*(x**2 + y**2)) - &
        14*mu*Rgas*u0**2*Sin(2*(x**2 + y**2)) + &
        12*k*v0**2*Sin(2*(x**2 + y**2)) - 12*g*k*v0**2*Sin(2*(x**2 + y**2)) + &
        14*mu*Rgas*v0**2*Sin(2*(x**2 + y**2)) + &
        18*et0*Rgas*rho0*u0*x*Sin(2*(x**2 + y**2)) - &
        16*mu*Rgas*u0*v0*x**2*Sin(2*(x**2 + y**2)) - &
        18*et0*Rgas*rho0*v0*y*Sin(2*(x**2 + y**2)) + &
        16*mu*Rgas*u0*v0*x*y*Sin(2*(x**2 + y**2)) - &
        16*mu*Rgas*u0*v0*y**2*Sin(2*(x**2 + y**2)) + &
        9*et0*Rgas*rho0*u0*x*Sin(3*(x**2 + y**2)))/(6.*Rgas)
    end function s_energy

end module mms
