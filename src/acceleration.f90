!===========================================================
! This module contains subroutines for accelerating the
! solution.
!===========================================================

module acceleration
  use mesh_class, only : mesh,element
  use linalg,     only : thomas
  implicit none
  private :: compute_mach_i,compute_mach_j
  public  :: irs_upwind

  contains

    !--------------------------------------------------------
    ! Function for smoothing the residual using central
    ! implicit residual smoothing
    !--------------------------------------------------------
    !pure function irs_central() result(r_bar)
    !end function irs_central

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
      p = (g-1.0d0)*(elem%u(4) - 0.5d0*rho*(u**2 + v**2))
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
      p = (g-1.0d0)*(elem%u(4) - 0.5d0*rho*(u**2 + v**2))
      a = sqrt(g*p/rho)

      mach = uproj/a

    end function compute_mach_j

    !----------------------------------------------------------------------------
    ! Function for smoothing the residual using upwind
    ! implicit residual smoothing
    !
    ! Notes:
    ! - The subroutine smooths first in the i-direction, then
    !   deallocates and reallocates memory for the
    !   the j-direction.  The i-smoothed residuals are then
    !   smoothed in the j-direction.
    !----------------------------------------------------------------------------
    subroutine irs_upwind(grid,eps,resid,g,r_bar2)
      implicit none
      type(mesh), intent(in)                            :: grid
      double precision, intent(in)                      :: eps
      double precision, intent(in), dimension(:,:,:)    :: resid
      double precision, intent(in)                      :: g
      double precision, intent(inout), dimension(:,:,:) :: r_bar2(:,:,:)
      double precision, allocatable                     :: r_bar1(:,:,:)
      double precision, allocatable                     :: mach(:,:)
      double precision, allocatable, dimension(:)       :: ld,d,ud,rhs,rslice
      double precision                                  :: eps_i,eps_j
      integer                                           :: i,j,k,aer,rshape(3)

      ! Allocating memory for the smoothed residual
      allocate(r_bar1,mold=resid,stat=aer)
      if (aer.ne.0) then
        print *, "Error: Can't allocate memory for r_bar1 in irs_upwind."
        stop
      end if

      ! Allocating memory for the Mach number
      rshape = shape(resid)
      allocate(mach(rshape(1),rshape(2)),stat=aer)
      if (aer.ne.0) then
        print *, "Error: Can't allocate memory for Mach array in irs_upwind."
        stop
      end if

      ! Allocating memory for inputs to thomas algorithm -- only for i-direction
      allocate(d(grid%nelemi))
      allocate(ld(grid%nelemi))
      allocate(ud(grid%nelemi-1))
      allocate(rhs(grid%nelemi))
      allocate(rslice(grid%nelemi))

      ! Computing the Mach number in the i-direction in each element
      do j=2,grid%nelemj-1
        do i=2,grid%nelemi-1
          mach(i,j) = compute_mach_i(grid%elem(i,j),g)
        end do
      end do

      ! Smoothing in the i-direction
      d(1)  = 1.0d0
      ud(1) = 0.0d0
      ld(grid%nelemi) = 0.0d0
      d(grid%nelemi)  = 1.0d0
      do k=1,4
        do j=2,grid%nelemj-1

          ! Setting the solution in the first and last elements
          rhs(1) = resid(1,j,k)
          rhs(grid%nelemi) = resid(grid%nelemi,j,k)

          ! Assembling system based on Mach number in i-direction
          do i=2,grid%nelemi-1

            ! Computing smoothing coefficient epsilon_i
            eps_i = eps*min(grid%elem(i,j)%lambda_ci/grid%elem(i,j)%lambda_cj,1.0d0)

            rhs(i) = resid(i,j,k)
            if (mach(i,j).gt.1.0d0) then
              ld(i) = -eps
              d(i)  = 1.0d0 + eps
              ud(i) = 0.0d0
            else if (abs(mach(i,j)).le.1.0d0) then
              ld(i) = -eps
              d(i)  = 1.0d0 + 2.0d0*eps
              ud(i) = -eps
            else  ! mach(i,j) < -1
              ld(i) = 0.0d0
              d(i)  = 1.0d0 + eps
              ud(i) = -eps
            end if
          end do

          ! Solving system using thomas algorithm
          call thomas(ud,ld,rhs,d,rslice)

          ! Copying solution into r_bar1
          r_bar1(:,j,k) = rslice

        end do
      end do

      ! De-allocating memory and reallocating for smoothing in j-direction
      deallocate(d)
      deallocate(ld)
      deallocate(ud)
      deallocate(rhs)
      deallocate(rslice)
      allocate(d(grid%nelemj))
      allocate(ld(grid%nelemj))
      allocate(ud(grid%nelemj-1))
      allocate(rhs(grid%nelemj))
      allocate(rslice(grid%nelemj))

      ! Computing the Mach number in the j-direction
      do j=2,grid%nelemj-1
        do i=2,grid%nelemi-1
          mach(i,j) = compute_mach_j(grid%elem(i,j),g)
        end do
      end do

      ! Smoothing in the j-direction
      d(1) = 1.0d0
      ud(1) = 0.0d0
      d(grid%nelemj) = 1.0d0
      ld(grid%nelemj) = 0.0d0
      do k=1,4
        do i=2,grid%nelemi-1

          ! Setting the solution in the first and last elements
          rhs(1) = r_bar1(i,1,k)
          rhs(grid%nelemj) = r_bar1(i,grid%nelemj,k)

          ! Assembling the system based on the Mach number in the j-direction
          do j=2,grid%nelemj-1

            ! Scaling the smoothing coefficient
            eps_j = eps*min(grid%elem(i,j)%lambda_cj/grid%elem(i,j)%lambda_ci,1.0d0)

            rhs(j) = r_bar1(i,j,k)
            if (mach(i,j).gt.1.0d0) then
              ld(j)  = -eps
              d(j)   = 1.0d0 + eps
              ud(j)  = 0.0d0
            else if (abs(mach(i,j)).le.1.0d0) then
              ld(j)  = -eps
              d(j)   = 1.0d0 + 2.0d0*eps
              ud(j)  = -eps
            else  ! mach(i,j) < -1
              ld(j)  = 0.0d0
              d(j)   = 1.0d0 + eps
              ud(j)  = -eps
            end if
          end do

          ! Solving system of equations
          call thomas(ud,ld,rhs,d,rslice)

          ! Copying the solution into r_bar2
          r_bar2(i,:,k) = rslice

        end do
      end do

    end subroutine irs_upwind

end module acceleration
