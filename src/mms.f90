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
  public :: rho_e,u_e,v_e,et_e,s_continuity,s_xmom,s_ymom,s_energy

  ! Parameters for MMS
  double precision, parameter :: rho0 = 0.5d0
  double precision, parameter :: u0   = 1.0d0
  double precision, parameter :: v0   = 0.1d0
  double precision, parameter :: et0  = 0.5d0
  double precision, parameter :: g    = 1.4d0
  double precision, parameter :: R    = 1.0d0
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
      s = (6*((-1 + g)*rho0*(et0 - eps*v0**2)*x + u0*v0*y)*Cos(x**2 + y**2)**2 - &
    3*(-1 + g)*rho0*v0**2*x*Cos(x**2 + y**2)**3 + &
    Sin(x**2 + y**2)*(-9*(-1 + g)*rho0*(et0 - eps*v0**2)*x - &
       4*mu*(v0 - 4*u0*x**2 + 2*u0*x*y - 3*u0*y**2) - &
       6*(-1 + g)*rho0*(et0 - eps*v0**2)*x*Sin(x**2 + y**2)) + &
    Cos(x**2 + y**2)*(6*eps*u0**2*x - &
       3*(-1 + g)*rho0*(-3*et0 + eps*((3 + eps)*u0**2 + eps*v0**2))*x + &
       4*mu*v0*(x - 2*y)*y + 2*u0*(-7*mu + 3*eps*v0*y) + &
       3*x*Sin(x**2 + y**2)*(2*u0**2 - &
          (-1 + g)*rho0*((3 + 4*eps)*u0**2 - 3*v0**2) - &
          (-1 + g)*rho0*(3*u0**2 - 2*v0**2)*Sin(x**2 + y**2))))/3.
    end function s_xmom

    !-------------------------------------------------------
    ! Function for computing the source term for the
    ! y-momentum equation
    !-------------------------------------------------------
    pure function s_ymom(x,y) result(s)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: s
      s = (6*(-1 + g)*rho0*(et0 - eps*v0**2)*y*Cos(x**2 + y**2)**2 - &
    3*(-1 + g)*rho0*v0**2*y*Cos(x**2 + y**2)**3 + &
    Sin(x**2 + y**2)*(-9*(-1 + g)*rho0*(et0 - eps*v0**2)*y - &
       6*eps*v0*(u0*x + v0*y) + 2*mu*(7*v0 + 2*u0*x*(-2*x + y)) - &
       6*(u0*v0*x + (-1 + g)*rho0*(et0 - eps*v0**2)*y)*Sin(x**2 + y**2)) + &
    Cos(x**2 + y**2)*(3*(-1 + g)*rho0* &
        (3*et0 - eps*((3 + eps)*u0**2 + eps*v0**2))*y + &
       4*mu*(u0 + v0*(3*x**2 - 2*x*y + 4*y**2)) + &
       3*y*Sin(x**2 + y**2)*(-2*v0**2 + &
          (-1 + g)*rho0*(-((3 + 4*eps)*u0**2) + 3*v0**2) - &
          (-1 + g)*rho0*(3*u0**2 - 2*v0**2)*Sin(x**2 + y**2))))/3.
    end function s_ymom

    !-------------------------------------------------------
    ! Function for computing the source term for the
    ! energy equation
    !-------------------------------------------------------
    pure function s_energy(x,y) result(s)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: s
      s = ((3*et0*(R*rho0*(9*u0*x + 6*eps*u0*x + v0*y + 6*eps*v0*y) + &
          8*(-1 + g)*k*(x**2 + y**2)) + &
       4*eps*(6*(-1 + g)*k*(u0**2 - v0**2*(x**2 + y**2)) + &
          mu*R*(-7*u0**2 + 2*u0*v0*(1 + x*y - 2*y**2) + &
             2*v0**2*(3*x**2 - 2*x*y + 4*y**2))))*Cos(x**2 + y**2) + &
    2*(12*(-1 + g)*k*(u0**2 - v0**2)*(x**2 + y**2) + &
       R*(u0*(4*mu*v0 + 3*(3 + 2*eps)*et0*rho0*x) - &
          4*mu*u0**2*(4*x**2 - 2*x*y + 3*y**2) + &
          v0*(3*(3 + 2*eps)*et0*rho0*y + 4*mu*v0*(3*x**2 - 2*x*y + 4*y**2))))* &
            Cos(2*(x**2 + y**2)) + 9*et0*R*rho0*v0*y*Cos(3*(x**2 + y**2)) - &
    24*et0*k*Sin(x**2 + y**2) + 24*et0*g*k*Sin(x**2 + y**2) - &
    8*eps*mu*R*u0*v0*Sin(x**2 + y**2) + 24*eps*k*v0**2*Sin(x**2 + y**2) - &
    24*eps*g*k*v0**2*Sin(x**2 + y**2) + 28*eps*mu*R*v0**2*Sin(x**2 + y**2) - &
    3*et0*R*rho0*u0*x*Sin(x**2 + y**2) - &
    18*eps*et0*R*rho0*u0*x*Sin(x**2 + y**2) + &
    24*eps*k*u0**2*x**2*Sin(x**2 + y**2) - &
    24*eps*g*k*u0**2*x**2*Sin(x**2 + y**2) + &
    32*eps*mu*R*u0**2*x**2*Sin(x**2 + y**2) - &
    16*eps*mu*R*u0*v0*x**2*Sin(x**2 + y**2) - &
    27*et0*R*rho0*v0*y*Sin(x**2 + y**2) - &
    18*eps*et0*R*rho0*v0*y*Sin(x**2 + y**2) - &
    16*eps*mu*R*u0**2*x*y*Sin(x**2 + y**2) + &
    8*eps*mu*R*u0*v0*x*y*Sin(x**2 + y**2) + &
    24*eps*k*u0**2*y**2*Sin(x**2 + y**2) - &
    24*eps*g*k*u0**2*y**2*Sin(x**2 + y**2) + &
    24*eps*mu*R*u0**2*y**2*Sin(x**2 + y**2) - &
    12*k*u0**2*Sin(2*(x**2 + y**2)) + 12*g*k*u0**2*Sin(2*(x**2 + y**2)) - &
    14*mu*R*u0**2*Sin(2*(x**2 + y**2)) + 12*k*v0**2*Sin(2*(x**2 + y**2)) - &
    12*g*k*v0**2*Sin(2*(x**2 + y**2)) + 14*mu*R*v0**2*Sin(2*(x**2 + y**2)) + &
    18*et0*R*rho0*u0*x*Sin(2*(x**2 + y**2)) - &
    16*mu*R*u0*v0*x**2*Sin(2*(x**2 + y**2)) - &
    18*et0*R*rho0*v0*y*Sin(2*(x**2 + y**2)) + &
    16*mu*R*u0*v0*x*y*Sin(2*(x**2 + y**2)) - &
    16*mu*R*u0*v0*y**2*Sin(2*(x**2 + y**2)) + &
    9*et0*R*rho0*u0*x*Sin(3*(x**2 + y**2)))/(6.*R)
    end function s_energy

end module mms
