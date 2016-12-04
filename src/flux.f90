!===========================================================
! This module contains functions for solving the Riemann 
! problem.  It also contains a function which is used to 
! compute the actual value of the flux when provided with 
! the vector of primitive variables.
!===========================================================
module flux
  use utils, only : w_to_u
  implicit none
  private
  public :: fluxax, fluxay

  contains

    !------------------------------------------------------
    ! Function for computing the x-flux term for the Euler
    ! equations (i.e., the advective flux).
    !
    ! parameters:
    ! - w: vector of primitive variables.
    ! - g: ratio of specific heats
    !------------------------------------------------------
    function fluxax(w,g) result(f)
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
    function fluxay(w,g) result(f)
      implicit none
      double precision, intent(in) :: w(4),g
      double precision             :: f(4),u(4)
      u = w_to_u(w,g)
      f(1) = w(1)*w(3)
      f(2) = w(1)*w(2)*w(3)
      f(3) = w(1)*w(3)**2 + w(4)
      f(4) = w(3)*(u(4)+w(4))
    end function fluxay

end module flux
