!===========================================================
! This module contains subroutines for accelerating the
! solution.
!===========================================================

module acceleration
  implicit none
  use mesh_class, only : mesh,element
  use linalg,     only : thomas
  private :: compute_mach_i,compute_mach_j
  public  :: irs_central,irs_upwind

  !--------------------------------------------------------
  ! Function for smoothing the residual using central
  ! implicit residual smoothing
  !--------------------------------------------------------
  pure function irs_central() result(r_bar)
  end function irs_central

  !--------------------------------------------------------
  ! Function for computing the Mach number in the i-dir
  !--------------------------------------------------------
  pure function compute_mach_i(elem,g) result(mach)
    implicit none
    type(element), intent(in)      :: elem
    double precision, intent(in)   :: g
    double precision               :: mach
    double precision, dimension(2) :: nhat_i
    double precision               :: u,v,a,uproj,p,rho
    integer                        :: aer

    ! Computing the i-averaged normal vector for the element
    nhat_i = 0.5d0*(elem%n(:,2) - elem%n(:,4))

   ! Computing the velocity projected into the i-direction
    rho = elem%u(1)
    u = elem%u(2)/rho
    v = elem%u(3)/rho
    uproj = u*nhat_i(1) + v*nhat_i(2)

    ! Computing the speed of sound
    p = (g-1.0d0)*(u(4) - 0.5d0*(u**2 + v**2))
    a = sqrt(g*p/rho)

    mach = uproj/a

  end function compute_mach_i

  !--------------------------------------------------------
  ! Function for computing the Mach number in the j-dir
  !--------------------------------------------------------
  pure function compute_mach_j(elem,g) result(mach)
    implicit none
    type(element), intent(in)      :: elem
    double precision, intent(in)   :: g
    double precision               :: mach
    double precision, dimension(2) :: nhat_j
    double precision               :: u,v,a,uproj,p,rho
    integer                        :: i,j,aer

    ! Computing the i-averaged normal vector for the element
    nhat_j = 0.5d0*(elem%n(:,3) - elem%n(:,1))

    ! Computing the velocity projected into the j-direction
    rho = elem%u(1)
    u = elem%u(2)/rho
    v = elem%u(3)/rho
    uproj = u*nhat_j(1) + v*nhat_j(2)

    ! Computing the speed of sound
    p = (g-1.0d0)*(u(4) - 0.5d0*(u**2 + v**2))
    a = sqrt(g*p/rho)

    mach = uproj/a

  end function compute_mach_i

  !--------------------------------------------------------
  ! Function for smoothing the residual using upwind
  ! implicit residual smoothing
  !--------------------------------------------------------
  subroutine irs_upwind(grid,eps,resid,r_bar)
    implicit none
    type(mesh), intent(in)                         :: grid
    double precision, intent(in)                   :: eps
    double precision, intent(in), dimension(:,:,:) :: resid(:,:,:)
    double precision, intent(inout), allocatable   :: r_bar(:,:,:)
    integer                                        :: i,j,aer

    ! Allocating memory for the smoothed residual
    allocate(r_bar,mold=resid,stat=aer)
    if (aer.ne.0) then
      print *, "Can't allocate memory for r_bar in irs_upwind."
      stop
    end if

    ! Computing the Mach number in the i-direction
    do j=1,grid%nelemj
      do i=1,grid%nelemi
      end do
    end do



  end function irs_upwind

end module acceleration
