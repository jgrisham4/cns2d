!===========================================================
! This module contains functions for solving the Riemann
! problem.  It also contains a function which is used to
! compute the actual value of the flux when provided with
! the vector of primitive variables.
!===========================================================
module flux
  use mesh_class,     only : element
  use utils,          only : w_to_u,u_to_w
  use gas_properties, only : mu,k
  implicit none
  private
  public :: fluxax, fluxay, flux_adv, flux_visc, flux_visc_state

  contains

    !------------------------------------------------------
    ! Function for computing the x-flux term for the Euler
    ! equations (i.e., the advective flux).
    !
    ! parameters:
    ! - w: vector of primitive variables.
    ! - g: ratio of specific heats
    !------------------------------------------------------
    pure function fluxax(w,g) result(f)
      implicit none
      double precision, intent(in) :: w(4),g
      double precision             :: f(4),u(4)
      u = w_to_u(w,g)
      f(1) = w(1)*w(2)
      f(2) = w(1)*w(2)**2 + w(4)
      f(3) = w(1)*w(2)*w(3)
      f(4) = w(2)*(u(4)+w(4))
    end function fluxax

    !------------------------------------------------------
    ! Function for computing the y-flux term for the Euler
    ! equations (i.e., the advective flux).
    !
    ! parameters:
    ! - w: vector of primitive variables.
    ! - g: ratio of specific heats
    !------------------------------------------------------
    pure function fluxay(w,g) result(f)
      implicit none
      double precision, intent(in) :: w(4),g
      double precision             :: f(4),u(4)
      u = w_to_u(w,g)
      f(1) = w(1)*w(3)
      f(2) = w(1)*w(2)*w(3)
      f(3) = w(1)*w(3)**2 + w(4)
      f(4) = w(3)*(u(4)+w(4))
    end function fluxay

    !------------------------------------------------------
    ! Function for computing the advective fluxes along
    ! a face.
    !------------------------------------------------------
    pure function flux_adv(w,n,g) result(f)
      implicit none
      double precision, intent(in) :: w(4),n(2),g
      double precision             :: f(4),u(4)
      double precision             :: rho,vx,vy,p,ht,vc

      ! Converting to conservative variables
      u = w_to_u(w,g)

      ! Separating variables
      rho = w(1)
      vx  = w(2)
      vy  = w(3)
      p   = w(4)
      ht  = u(4) + p/rho  ! ht = et + p/rho

      ! Finding contravariant velocity
      vc = vx*n(1) + vy*n(2)

      ! Forming flux
      f(1) = rho*vc
      f(2) = rho*vx*vc + n(1)*p
      f(3) = rho*vy*vc + n(2)*p
      f(4) = rho*ht*vc

    end function flux_adv

    !------------------------------------------------------
    ! Function for computing the viscous fluxes along
    ! a face when the state at the face is provided.
    !------------------------------------------------------
    pure function flux_visc_state(u,v,T,dudx,dudy,dvdx,dvdy,dTdx,dTdy,n,g,R) result(f)
      implicit none
      double precision, intent(in) :: u,v,T,dudx,dudy,dvdx,dvdy,dTdx,dTdy,n(2)
      double precision, intent(in) :: g,R
      double precision             :: f(4),m,kval,lambda
      double precision             :: tauxx,tauxy,tauyy
      double precision             :: thetax,thetay

      ! Using Stokes' Theorem to compute lambda
      m     = mu(T)
      kval  = k(T,g,R)
      lambda = -2.0d0/3.0d0*m

      ! Computing components of shear stress tensor
      tauxx = lambda*(dudx + dvdy) + 2.0d0*m*dudx
      tauyy = lambda*(dudx + dvdy) + 2.0d0*m*dvdy
      tauxy = m*(dudy + dvdx)

      ! Computing work done by shear in addition to heat transfer
      thetax = u*tauxx + v*tauxy + kval*dTdx
      thetay = u*tauxy + v*tauyy + kval*dTdy

      ! Computing the flux
      f(1) = 0.0d0
      f(2) = n(1)*tauxx  + n(2)*tauxy
      f(3) = n(1)*tauxy  + n(2)*tauyy
      f(4) = n(1)*thetax + n(2)*thetay

    end function flux_visc_state

    !------------------------------------------------------
    ! Function for computing the viscous fluxes along
    ! a face when the elements on either side of the
    ! interface are provided as inputs.
    !------------------------------------------------------
    pure function flux_visc(elemL,elemR,n,g,R) result(f)
      implicit none
      type(element), intent(in)    :: elemL,elemR
      double precision, intent(in) :: n(2),g,R
      double precision             :: f(4),rho
      double precision             :: u,v,T,TL,TR,wL(4),wR(4)
      double precision             :: dudx,dudy,dvdx,dvdy,dTdx,dTdy

      ! Converting from conservative variables to primitive
      wL = u_to_w(elemL%u,g)
      wR = u_to_w(elemR%u,g)

      ! Computing temperature
      TL = wL(4)/(wL(1)*R)
      TR = wR(4)/(wR(1)*R)

      ! Computing averaged properties at the interface
      u    = 0.5d0*(wL(2) + wR(2))
      v    = 0.5d0*(wL(3) + wR(3))
      T    = 0.5d0*(TL + TR)
      dudx = 0.5d0*(elemL%dVdx(1) + elemR%dVdx(1))
      dudy = 0.5d0*(elemL%dVdy(1) + elemR%dVdy(1))
      dvdx = 0.5d0*(elemL%dVdx(2) + elemR%dVdx(2))
      dvdy = 0.5d0*(elemL%dVdy(2) + elemR%dVdy(2))
      dTdx = 0.5d0*(elemL%dTdx    + elemR%dTdx)
      dTdy = 0.5d0*(elemL%dTdy    + elemR%dTdy)

      ! Calling function to compute the viscous fluxes
      f = flux_visc_state(u,v,T,dudx,dudy,dvdx,dvdy,dTdx,dTdy,n,g,R)

    end function flux_visc

end module flux
