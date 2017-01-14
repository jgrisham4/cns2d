!===============================================================================
! This module contains the subroutines for a simple Euler solver.
!
! Author: James Grisham
! Date: 01/13/2017
!===============================================================================

module euler
  use solver_class, only : solver
  use utils,        only : compute_elem_max,compute_elem_min
  use limiters,     only : barth
  use grad,         only : compute_gradient
  use riemann,      only : roe,rotated_rhll
  use bcs,          only : apply_bcs
  implicit none
  public :: residual_inv,residual_inv_fo

  contains

    !---------------------------------------------------------------------------
    ! Subroutine for computing the residual for Euler eqs
    !---------------------------------------------------------------------------
    subroutine residual_inv(this,resid)
      implicit none
      type(solver), intent(inout)     :: this
      double precision, intent(inout) :: resid(:,:,:)
      double precision                :: duL(4),duR(4),phi(4)
      double precision, dimension(4)  :: fx,fy,uextrap,wextrap,wtmp
      double precision, dimension(2)  :: rL,rR,r
      double precision, allocatable   :: u(:,:,:),gradU(:,:,:),umax(:,:,:),umin(:,:,:)
      integer                         :: i,j,k,l,err

      ! Allocating memory for gradient of state at cell centers
      allocate(gradU(this%grid%nelemi,this%grid%nelemj,8),stat=err)
      if (err.ne.0) then
        print *, "Error: Can't allocate memory for gradU in residual_inv."
        stop
      end if
      allocate(u(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for u in residual_inv."
        stop
      end if
      allocate(umax(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for umax in residual_inv."
        stop
      end if
      allocate(umin(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for umin in residual_inv."
        stop
      end if

      ! Copying the last solution into u
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          u(i,j,:) = this%grid%elem(i,j)%u
        end do
      end do

      ! Computing gradient for use in reconstruction
      call compute_gradient(this%grid,u,gradU)

      ! Assigning gradient to elements
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%elem(i,j)%dudx = gradU(i,j,1:4)
          this%grid%elem(i,j)%dudy = gradU(i,j,5:8)
        end do
      end do

      ! Finding max and min of states for each element and neighbors
      if (this%limiter.eq."barth") then
        umax = compute_elem_max(u,this%grid%nelemi,this%grid%nelemj)
        umin = compute_elem_min(u,this%grid%nelemi,this%grid%nelemj)
      end if

      ! Reconstructing states on left and right sides of each vertical interface
      ! for each element
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi

          if (i.eq.1) then

            ! Only need to reconstruct the left state
            rL(1) = this%grid%edges_v(i+1,j)%xm - this%grid%elem(i,j)%xc
            rL(2) = this%grid%edges_v(i+1,j)%ym - this%grid%elem(i,j)%yc

            ! Finding grad(u) . r_L on the left side of the interface
            duL = gradU(i,j,1:4)*rL(1) + gradU(i,j,5:8)*rL(2)

            ! Reconstructing state on left of interface
            select case (this%limiter)
              case ("none")
                do k=1,4
                  this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)
                end do
              case ("barth")

                ! Calling barth subroutine to find limiter
                phi = barth(this,i,j,gradU,umax,umin)

                ! Slope-limited reconstruction
                do k=1,4
                  this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)*phi(k)
                end do

              case default
                write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
                stop
            end select

          else if (i.eq.(this%grid%nelemi)) then

            ! Only need to reconstruct the right state
            rR(1) = this%grid%edges_v(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_v(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding the slope on the right side of the interface
            duR = gradU(i,j,1:4)*rR(1) + gradU(i,j,5:8)*rR(2)

            ! Reconstructing primitive states on right of interface
            select case (this%limiter)
              case ("none")
                do k=1,4
                  this%grid%edges_v(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + duR(k)
                end do
              case ("barth")

                ! Calling subroutine to compute the barth limiter
                phi = barth(this,i,j,gradU,umax,umin)

                ! Slope-limited reconstruction
                do k=1,4
                  this%grid%edges_v(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + duR(k)*phi(k)
                end do

              case default
                write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
                stop
            end select

          else

            ! Computing position vectors for face midpoints
            rL(1) = this%grid%edges_v(i+1,j)%xm - this%grid%elem(i,j)%xc
            rL(2) = this%grid%edges_v(i+1,j)%ym - this%grid%elem(i,j)%yc
            rR(1) = this%grid%edges_v(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_v(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding slopes on left and right
            ! The below is grad(u) . r_L and grad(u) . r_R
            duL = gradU(i,j,1:4)*rL(1) + gradU(i,j,5:8)*rL(2)
            duR = gradU(i,j,1:4)*rR(1) + gradU(i,j,5:8)*rR(2)

            ! Reconstructing primitive states on left and right of interface
            select case (this%limiter)
              case ("none")
                do k=1,4
                  this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)
                  this%grid%edges_v(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + duR(k)
                end do
              case ("barth")

                ! Calling subroutine to compute the barth limiter
                phi = barth(this,i,j,gradU,umax,umin)

                ! Slope-limited reconstruction
                do k=1,4
                  this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)*phi(k)
                  this%grid%edges_v(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + duR(k)*phi(k)
                end do

              case default
                write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
                stop
            end select

          end if

        end do
      end do

      ! Reconstructing states on left and right sides of each horizontal interface
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi

          if (j.eq.1) then

            ! Only need to reconstruct the left state
            rL(1) = this%grid%edges_h(i,j+1)%xm - this%grid%elem(i,j)%xc
            rL(2) = this%grid%edges_h(i,j+1)%ym - this%grid%elem(i,j)%yc

            ! Finding the slope on the left side of the interface
            duL = gradU(i,j,1:4)*rL(1) + gradU(i,j,5:8)*rL(2)

            ! Reconstructing primitive states on left of interface
            select case (this%limiter)
              case ("none")
                do k=1,4
                  this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)
                end do
              case ("barth")

                ! Calling barth subroutine to find limiter
                phi = barth(this,i,j,gradU,umax,umin)

                ! Slope-limited reconstruction
                do k=1,4
                  this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)*phi(k)
                end do

              case default
                write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
                stop
            end select

          else if (j.eq.(this%grid%nelemj)) then

            ! Only need to reconstruct the right state
            rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding the slope on the right side of the interface
            duR = gradU(i,j,1:4)*rR(1) + gradU(i,j,5:8)*rR(2)

            ! Reconstructing primitive states on right of interface
            select case (this%limiter)
              case ("none")
                do k=1,4
                  this%grid%edges_h(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + duR(k)
                end do
              case ("barth")

                ! Calling subroutine to compute the barth limiter
                phi = barth(this,i,j,gradU,umax,umin)

                ! Slope-limited reconstruction
                do k=1,4
                  this%grid%edges_h(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + duR(k)*phi(k)
                end do

              case default
                write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
                stop
            end select

          else

            ! Computing position vectors for face midpoints
            rL(1) = this%grid%edges_h(i,j+1)%xm - this%grid%elem(i,j)%xc
            rL(2) = this%grid%edges_h(i,j+1)%ym - this%grid%elem(i,j)%yc
            rR(1) = this%grid%edges_h(i,j)%xm   - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym   - this%grid%elem(i,j)%yc

            ! Finding slopes on left and right
            duL = gradU(i,j,1:4)*rL(1) + gradU(i,j,5:8)*rL(2)
            duR = gradU(i,j,1:4)*rR(1) + gradU(i,j,5:8)*rR(2)

            ! Reconstructing primitive states on left and right of interface
            select case (this%limiter)
              case ("none")
                do k=1,4
                  this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)
                  this%grid%edges_h(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + duR(k)
                end do
              case ("barth")

                ! Calling subroutine to compute the barth limiter
                phi = barth(this,i,j,gradU,umax,umin)

                ! Slope-limited reconstruction
                do k=1,4
                  this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)*phi(k)
                  this%grid%edges_h(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + duR(k)*phi(k)
                end do

              case default
                write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
                stop
            end select

          end if

        end do
      end do

      ! Applying boundary conditions
      call apply_bcs(this)

      ! Must now iterate through all the interior interfaces and solve
      ! the Riemann problem to find the fluxes
      ! Vertical faces
      do j=1,this%grid%nelemj
        do i=2,this%grid%nelemi

          ! Calling Riemann solver
          !print *, i, j, this%grid%edges_v(i,j)%xm, this%grid%edges_v(i,j)%ym
          this%grid%edges_v(i,j)%flux = roe(this%grid%edges_v(i,j)%uL,this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2))

          !this%grid%edges_v(i,j)%flux = rotated_rhll(this%grid%edges_v(i,j)%uL,this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2))

          ! Checking for NaNs
          !do k=1,4
          !  if (isnan(this%grid%edges_v(i,j)%flux(k))) then
          !    write(*,'(a)') "NaNs encountered after solving Riemann problem for vertical faces"
          !    write(*,'(a,i4,a,i4,a,i4)') "i = ", i, " j = ", j , " k = ", k
          !    write(*,'(2(a,f12.5))') "x_midpoint = ", this%grid%edges_v(i,j)%xm, &
          !      " y_midpoint = ", this%grid%edges_v(i,j)%ym
          !    stop
          !  end if
          !end do

        end do
      end do

      ! Horizontal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi

          ! Calling Riemann solver
          this%grid%edges_h(i,j)%flux = roe(this%grid%edges_h(i,j)%uL,this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3))
          !this%grid%edges_h(i,j)%flux = rotated_rhll(this%grid%edges_h(i,j)%uL,this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3))

          ! Checking for NaNs
          !do k=1,4
          !  if (isnan(this%grid%edges_h(i,j)%flux(k))) then
          !    write(*,'(a)') "NaNs encountered after solving Riemann problem for horizontal faces"
          !    write(*,'(a,i4,a,i4,a,i4)') "i = ", i, " j = ", j , " k = ", k
          !    write(*,'(2(a,f12.5))') "x_midpoint = ", this%grid%edges_h(i,j)%xm, &
          !      " y_midpoint = ", this%grid%edges_h(i,j)%ym
          !    stop
          !  end if
          !end do

        end do
      end do

      ! Computing residual
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          resid(i,j,:) = -1.0d0/(this%grid%elem(i,j)%area)* &
            (-this%grid%edges_h(i,j)%flux*this%grid%edges_h(i,j)%length + &
            this%grid%edges_h(i,j+1)%flux*this%grid%edges_h(i,j+1)%length - &
            this%grid%edges_v(i,j)%flux*this%grid%edges_v(i,j)%length + &
            this%grid%edges_v(i+1,j)%flux*this%grid%edges_v(i+1,j)%length)
        end do
      end do

    end subroutine residual_inv

    !---------------------------------------------------------------------------
    ! Subroutine for computing the residual for Euler eqs
    ! This subroutine is modified from the above to compute
    ! the first-order accurate residual.  This is accomplished
    ! by removing the piecewise linear reconstruction at the
    ! interface and replacing it with constant extrapolation
    ! or zero-th order interpolation.
    !---------------------------------------------------------------------------
    subroutine residual_inv_fo(this,resid)
      implicit none
      type(solver), intent(inout)     :: this
      double precision, intent(inout) :: resid(:,:,:)
      double precision                :: duL(4),duR(4),phi(4)
      double precision, dimension(4)  :: fx,fy,uextrap,wextrap,wtmp
      double precision, dimension(2)  :: rL,rR,r
      double precision, allocatable   :: u(:,:,:),gradU(:,:,:),umax(:,:,:),umin(:,:,:)
      integer                         :: i,j,k,l,err

      ! Allocating memory for gradient of state at cell centers
      allocate(gradU(this%grid%nelemi,this%grid%nelemj,8),stat=err)
      if (err.ne.0) then
        print *, "Error: Can't allocate memory for gradU in residual_inv_fo."
        stop
      end if
      allocate(u(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for u in residual_inv_fo."
        stop
      end if
      allocate(umax(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for umax in residual_inv_fo."
        stop
      end if
      allocate(umin(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for umin in residual_inv_fo."
        stop
      end if

      ! Copying the last solution into u
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          u(i,j,:) = this%grid%elem(i,j)%u
        end do
      end do

      ! Applying boundary conditions
      call apply_bcs(this)

      ! Must now iterate through all the interior interfaces and solve
      ! the Riemann problem to find the fluxes
      ! Vertical faces
      do j=1,this%grid%nelemj
        do i=2,this%grid%nelemi

          ! Calling Riemann solver
          !print *, i, j, this%grid%edges_v(i,j)%xm, this%grid%edges_v(i,j)%ym
          this%grid%edges_v(i,j)%flux = roe(this%grid%elem(i-1,j)%u,this%grid%elem(i,j)%u,this%grid%elem(i-1,j)%n(:,2))
          !this%grid%edges_v(i,j)%flux = rotated_rhll(this%grid%edges_v(i,j)%uL,this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2))

        end do
      end do

      ! Horizontal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi

          ! Calling Riemann solver
          this%grid%edges_h(i,j)%flux = roe(this%grid%elem(i,j-1)%u,this%grid%elem(i,j)%u,this%grid%elem(i,j-1)%n(:,3))
          !this%grid%edges_h(i,j)%flux = rotated_rhll(this%grid%edges_h(i,j)%uL,this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3))

        end do
      end do

      ! Computing residual
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          resid(i,j,:) = -1.0d0/(this%grid%elem(i,j)%area)* &
            (-this%grid%edges_h(i,j)%flux*this%grid%edges_h(i,j)%length + &
            this%grid%edges_h(i,j+1)%flux*this%grid%edges_h(i,j+1)%length - &
            this%grid%edges_v(i,j)%flux*this%grid%edges_v(i,j)%length + &
            this%grid%edges_v(i+1,j)%flux*this%grid%edges_v(i+1,j)%length)
        end do
      end do

    end subroutine residual_inv_fo

end module euler
