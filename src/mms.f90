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

  ! Defining pi
  double precision, parameter :: pi = 4.0d0*atan(1.0d0)

  ! Parameters for MMS
  double precision, parameter :: r0   = 1.0d0
  double precision, parameter :: u0   = 800.0d0
  double precision, parameter :: v0   = 800.0d0
  double precision, parameter :: p0   = 100000.0d0
  double precision, parameter :: rx   = 0.15d0
  double precision, parameter :: ux   = 50.0d0
  double precision, parameter :: vx   = -75.0d0
  double precision, parameter :: px   = 2000.0d0
  double precision, parameter :: ry   = -0.1d0
  double precision, parameter :: uy   = -30.0d0
  double precision, parameter :: vy   = 40.0d0
  double precision, parameter :: py   = 5000.0d0
  double precision, parameter :: arx  = 1.0d0
  double precision, parameter :: aux  = 1.5d0
  double precision, parameter :: avx  = 0.5d0
  double precision, parameter :: apx  = 2.0d0/3.0d0
  double precision, parameter :: ary  = 0.5d0
  double precision, parameter :: auy  = 0.6d0
  double precision, parameter :: avy  = 1.5d0
  double precision, parameter :: apy  = 1.0d0

  ! Gas properties
  double precision, parameter :: g    = 1.4d0
  double precision, parameter :: Rgas = 287.0d0
  double precision, parameter :: L    = 1.0d0
  double precision, parameter :: Pr   = 1.0d0
  double precision, parameter :: mu   = 20.0d0
  double precision, parameter :: k    = (mu*g*Rgas/(g-1.0d0))/Pr

  contains

    !-------------------------------------------------------
    ! Function for rho(x,y)
    !-------------------------------------------------------
    pure function rho_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L)
    end function rho_e

    !-------------------------------------------------------
    ! Function for u(x,y)
    !-------------------------------------------------------
    pure function u_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L)
    end function u_e

    !-------------------------------------------------------
    ! Function for v(x,y)
    !-------------------------------------------------------
    pure function v_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L)
    end function v_e

    !-------------------------------------------------------
    ! Function for et(x,y)
    !-------------------------------------------------------
    pure function et_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = (Rgas*(p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L)))/(-1 + g) + &
        ((u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L))**2 + &
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L))**2)/2.
    end function et_e

    !-------------------------------------------------------
    ! Function for computing du/dx
    !-------------------------------------------------------
    pure function dudx_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = (aux*pi*ux*Cos((aux*pi*x)/L))/L
    end function dudx_e

    !-------------------------------------------------------
    ! Function for computing du/dy
    !-------------------------------------------------------
    pure function dudy_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = -((auy*pi*uy*Sin((auy*pi*y)/L))/L)
    end function dudy_e

    !-------------------------------------------------------
    ! Function for computing dv/dx
    !-------------------------------------------------------
    pure function dvdx_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = -((avx*pi*vx*Sin((avx*pi*x)/L))/L)
    end function dvdx_e

    !-------------------------------------------------------
    ! Function for computing dv/dy
    !-------------------------------------------------------
    pure function dvdy_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = (avy*pi*vy*Cos((avy*pi*y)/L))/L
    end function dvdy_e

    !-------------------------------------------------------
    ! Function for computing dT/dx
    !-------------------------------------------------------
    pure function dTdx_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = -((apx*pi*px*Sin((apx*pi*x)/L))/&
        (L*Rgas*(r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L)))) - &
        (arx*pi*rx*Cos((arx*pi*x)/L)*&
        (p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L)))/&
        (L*Rgas*(r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))**2)
    end function dTdx_e

    !-------------------------------------------------------
    ! Function for computing dT/dy
    !-------------------------------------------------------
    pure function dTdy_e(x,y) result(r)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: r
      r = (apy*pi*py*Cos((apy*pi*y)/L))/&
        (L*Rgas*(r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))) + &
        (ary*pi*ry*(p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L))*&
        Sin((ary*pi*y)/L))/&
        (L*Rgas*(r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))**2)
    end function dTdy_e

    !-------------------------------------------------------
    ! Function for computing the source term for the
    ! continuity equation
    !-------------------------------------------------------
    pure function s_continuity(x,y) result(s)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: s
      s = (pi*(aux*ux*Cos((aux*pi*x)/L)*(r0 + ry*Cos((ary*pi*y)/L) + &
        rx*Sin((arx*pi*x)/L)) + &
        avy*vy*Cos((avy*pi*y)/L)*&
        (r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L)) + &
        arx*rx*Cos((arx*pi*x)/L)*&
        (u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L)) - &
        ary*ry*Sin((ary*pi*y)/L)*&
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L))))/L
    end function s_continuity

    !-------------------------------------------------------
    ! Function for computing the source term for the
    ! x-momentum equation
    !-------------------------------------------------------
    pure function s_xmom(x,y) result(s)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: s

      s = (auy**2*mu*Pi**2*uy*Cos((auy*Pi*y)/L) - apx*L*Pi*px*Sin((apx*Pi*x)/L) - &
        6.579736267392905*aux**2*mu*ux*Sin((aux*Pi*x)/L) + &
        2*aux**2*mu*Pi**2*ux*Sin((aux*Pi*x)/L) + &
        2*aux*L*Pi*ux*Cos((aux*Pi*x)/L)*&
        (r0 + ry*Cos((ary*Pi*y)/L) + rx*Sin((arx*Pi*x)/L))*&
        (u0 + uy*Cos((auy*Pi*y)/L) + ux*Sin((aux*Pi*x)/L)) + &
        avy*L*Pi*vy*Cos((avy*Pi*y)/L)*&
        (r0 + ry*Cos((ary*Pi*y)/L) + rx*Sin((arx*Pi*x)/L))*&
        (u0 + uy*Cos((auy*Pi*y)/L) + ux*Sin((aux*Pi*x)/L)) + &
        arx*L*Pi*rx*Cos((arx*Pi*x)/L)*&
        (u0 + uy*Cos((auy*Pi*y)/L) + ux*Sin((aux*Pi*x)/L))**2 - &
        ary*L*Pi*ry*(u0 + uy*Cos((auy*Pi*y)/L) + ux*Sin((aux*Pi*x)/L))*&
        Sin((ary*Pi*y)/L)*(v0 + vx*Cos((avx*Pi*x)/L) + vy*Sin((avy*Pi*y)/L)) - &
        auy*L*Pi*uy*(r0 + ry*Cos((ary*Pi*y)/L) + rx*Sin((arx*Pi*x)/L))*&
        Sin((auy*Pi*y)/L)*(v0 + vx*Cos((avx*Pi*x)/L) + vy*Sin((avy*Pi*y)/L)))/L**2

    end function s_xmom

    !-------------------------------------------------------
    ! Function for computing the source term for the
    ! y-momentum equation
    !-------------------------------------------------------
    pure function s_ymom(x,y) result(s)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: s

      s = (avx**2*mu*Pi**2*vx*Cos((avx*Pi*x)/L) + apy*L*Pi*py*Cos((apy*Pi*y)/L) - &
        avx*L*Pi*vx*(r0 + ry*Cos((ary*Pi*y)/L) + rx*Sin((arx*Pi*x)/L))*&
        (u0 + uy*Cos((auy*Pi*y)/L) + ux*Sin((aux*Pi*x)/L))*Sin((avx*Pi*x)/L) - &
        6.579736267392905*avy**2*mu*vy*Sin((avy*Pi*y)/L) + &
        2*avy**2*mu*Pi**2*vy*Sin((avy*Pi*y)/L) + &
        aux*L*Pi*ux*Cos((aux*Pi*x)/L)*&
        (r0 + ry*Cos((ary*Pi*y)/L) + rx*Sin((arx*Pi*x)/L))*&
        (v0 + vx*Cos((avx*Pi*x)/L) + vy*Sin((avy*Pi*y)/L)) + &
        2*avy*L*Pi*vy*Cos((avy*Pi*y)/L)*&
        (r0 + ry*Cos((ary*Pi*y)/L) + rx*Sin((arx*Pi*x)/L))*&
        (v0 + vx*Cos((avx*Pi*x)/L) + vy*Sin((avy*Pi*y)/L)) + &
        arx*L*Pi*rx*Cos((arx*Pi*x)/L)*&
        (u0 + uy*Cos((auy*Pi*y)/L) + ux*Sin((aux*Pi*x)/L))*&
        (v0 + vx*Cos((avx*Pi*x)/L) + vy*Sin((avy*Pi*y)/L)) - &
        ary*L*Pi*ry*Sin((ary*Pi*y)/L)*&
        (v0 + vx*Cos((avx*Pi*x)/L) + vy*Sin((avy*Pi*y)/L))**2)/L**2

    end function s_ymom

    !-------------------------------------------------------
    ! Function for computing the source term for the
    ! energy equation
    !-------------------------------------------------------
    pure function s_energy(x,y) result(s)
      implicit none
      double precision, intent(in) :: x,y
      double precision :: s

      s = (-13.159472534785811*aux*mu*ux*Cos((aux*pi*x)/L)*&
        (1.*aux*ux*Cos((aux*pi*x)/L) - &
        0.49999999999999994*avy*vy*Cos((avy*pi*y)/L)) - &
        13.159472534785811*avy*mu*vy*Cos((avy*pi*y)/L)*&
        (-0.49999999999999994*aux*ux*Cos((aux*pi*x)/L) + &
        1.*avy*vy*Cos((avy*pi*y)/L)) + &
        auy**2*mu*pi**2*uy*Cos((auy*pi*y)/L)*&
        (u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L)) - &
        apx*L*pi*px*Sin((apx*pi*x)/L)*&
        (u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L)) + &
        13.159472534785811*aux**2*mu*ux*Sin((aux*pi*x)/L)*&
        (u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L)) + &
        aux*L*pi*ux*Cos((aux*pi*x)/L)*&
        (p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L)) + &
        avy*L*pi*vy*Cos((avy*pi*y)/L)*&
        (p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L)) - &
        (k*pi**2*(2*apx*arx*px*rx*Cos((arx*pi*x)/L)*Sin((apx*pi*x)/L)*&
        (r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L)) - &
        apx**2*px*Cos((apx*pi*x)/L)*&
        (r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))**2 + &
        2*arx**2*rx**2*Cos((arx*pi*x)/L)**2*&
        (p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L)) + &
        arx**2*rx*Sin((arx*pi*x)/L)*&
        (r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))*&
        (p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L))))/&
        (Rgas*(r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))**3) - &
        (k*pi**2*(-(apy**2*py*(r0 + ry*Cos((ary*pi*y)/L) + &
        rx*Sin((arx*pi*x)/L))**2*Sin((apy*pi*y)/L)) + &
        ary**2*ry*Cos((ary*pi*y)/L)*&
        (r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))*&
        (p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L)) + &
        2*apy*ary*py*ry*Cos((apy*pi*y)/L)*&
        (r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))*&
        Sin((ary*pi*y)/L) + 2*ary**2*ry**2*&
        (p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L))*&
        Sin((ary*pi*y)/L)**2))/&
        (Rgas*(r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))**3) - &
        avx*mu*pi**2*vx*Sin((avx*pi*x)/L)*&
        (avx*vx*Sin((avx*pi*x)/L) + auy*uy*Sin((auy*pi*y)/L)) - &
        auy*mu*pi**2*uy*Sin((auy*pi*y)/L)*&
        (avx*vx*Sin((avx*pi*x)/L) + auy*uy*Sin((auy*pi*y)/L)) + &
        avx**2*mu*pi**2*vx*Cos((avx*pi*x)/L)*&
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L)) + &
        apy*L*pi*py*Cos((apy*pi*y)/L)*&
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L)) + &
        13.159472534785811*avy**2*mu*vy*Sin((avy*pi*y)/L)*&
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L)) + &
        L*pi*(r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))*&
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L))*&
        ((apy*py*Rgas*Cos((apy*pi*y)/L))/(-1 + g) - &
        auy*uy*(u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L))*&
        Sin((auy*pi*y)/L) + avy*vy*Cos((avy*pi*y)/L)*&
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L))) + &
        L*pi*(r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))*&
        (u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L))*&
        (-((apx*px*Rgas*Sin((apx*pi*x)/L))/(-1 + g)) + &
        aux*ux*Cos((aux*pi*x)/L)*&
        (u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L)) - &
        avx*vx*Sin((avx*pi*x)/L)*&
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L))) + &
        aux*L*pi*ux*Cos((aux*pi*x)/L)*&
        (r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))*&
        ((Rgas*(p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L)))/(-1 + g) + &
        ((u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L))**2 + &
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L))**2)/2.) + &
        avy*L*pi*vy*Cos((avy*pi*y)/L)*&
        (r0 + ry*Cos((ary*pi*y)/L) + rx*Sin((arx*pi*x)/L))*&
        ((Rgas*(p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L)))/(-1 + g) + &
        ((u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L))**2 + &
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L))**2)/2.) + &
        arx*L*pi*rx*Cos((arx*pi*x)/L)*&
        (u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L))*&
        ((Rgas*(p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L)))/(-1 + g) + &
        ((u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L))**2 + &
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L))**2)/2.) - &
        ary*L*pi*ry*Sin((ary*pi*y)/L)*&
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L))*&
        ((Rgas*(p0 + px*Cos((apx*pi*x)/L) + py*Sin((apy*pi*y)/L)))/(-1 + g) + &
        ((u0 + uy*Cos((auy*pi*y)/L) + ux*Sin((aux*pi*x)/L))**2 + &
        (v0 + vx*Cos((avx*pi*x)/L) + vy*Sin((avy*pi*y)/L))**2)/2.))/L**2

    end function s_energy

end module mms
