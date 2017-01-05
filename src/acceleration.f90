!===========================================================
! This module contains subroutines for accelerating the
! solution.
!===========================================================

module acceleration
  implicit none
  private
  public :: irs_central,irs_upwind

  !--------------------------------------------------------
  ! Function for smoothing the residual using central
  ! implicit residual smoothing
  !--------------------------------------------------------
  pure function irs_central() result(r_bar)
  end function irs_central

end module acceleration
