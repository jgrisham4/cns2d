!===============================================================================
! This module contains the subroutines for a simple Euler solver.
!
! Author: James Grisham
! Date: 01/13/2017
!===============================================================================

module euler
  use solver_class, only : solver
  use utils,        only : compute_elem_max,compute_elem_min
  use limiters,     only : barth,venkatakrishnan
  use grad,         only : compute_gradient
  use riemann,      only : roe,rotated_rhll
  use bcs,          only : apply_bcs
  implicit none
  public :: residual_inv,residual_inv_fo

  contains

    !---------------------------------------------------------------------------
    ! Subroutine for computing the residual for Euler eqs
    ! This subroutine only computes and returns the residual.  It does not
    ! update the solution.
    !---------------------------------------------------------------------------
    subroutine residual_inv(this,resid)
      implicit none
      type(solver),     intent(inout) :: this
      double precision, intent(inout) :: resid(:,:,:)
      double precision                :: duL(4),duR(4)
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

      ! Finding max and min of states for each element and neighbors
      if ((this%limiter=="barth").or.(this%limiter=="venkat")) then
        umax = compute_elem_max(u,this%grid%nelemi,this%grid%nelemj)
        umin = compute_elem_min(u,this%grid%nelemi,this%grid%nelemj)
      end if

      ! Computing value of limiter for each element and assigning gradients
      if (this%limiter=="barth") then
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%dudx = gradU(i,j,1:4)
            this%grid%elem(i,j)%dudy = gradU(i,j,5:8)
            this%grid%elem(i,j)%phi = barth(this,i,j,gradU,umax,umin)
          end do
        end do
      else if (this%limiter=="venkat") then
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%dudx = gradU(i,j,1:4)
            this%grid%elem(i,j)%dudy = gradU(i,j,5:8)
            this%grid%elem(i,j)%phi = venkatakrishnan(this,i,j,gradU,umax,umin)
          end do
        end do
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
            duL = this%grid%elem(i,j)%dudx*rL(1) + this%grid%elem(i,j)%dudy*rL(2)

            ! Reconstructing state on left of interface
            if (this%limiter=="none") then

              ! Unlimited reconstruction
              do k=1,4
                this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)
              end do

            else if ((this%limiter=="barth").or.(this%limiter=="venkat")) then

              ! Slope-limited reconstruction
              do k=1,4
                this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + &
                  duL(k)*this%grid%elem(i,j)%phi(k)
              end do

            else
              write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
              stop
            end if

          else if (i.eq.(this%grid%nelemi)) then

            ! Only need to reconstruct the right state
            rR(1) = this%grid%edges_v(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_v(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding the slope on the right side of the interface
            duR = this%grid%elem(i,j)%dudx*rR(1) + this%grid%elem(i,j)%dudy*rR(2)

            ! Reconstructing primitive states on right of interface
            if (this%limiter=="none") then

              ! Unlimited reconstruction
              do k=1,4
                this%grid%edges_v(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + duR(k)
              end do

            else if ((this%limiter=="barth").or.(this%limiter=="venkat")) then

              ! Slope-limited reconstruction
              do k=1,4
                this%grid%edges_v(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + &
                  duR(k)*this%grid%elem(i,j)%phi(k)
              end do

            else
              write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
              stop
            end if

          else

            ! Computing position vectors for face midpoints
            rL(1) = this%grid%edges_v(i+1,j)%xm - this%grid%elem(i,j)%xc
            rL(2) = this%grid%edges_v(i+1,j)%ym - this%grid%elem(i,j)%yc
            rR(1) = this%grid%edges_v(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_v(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding slopes on left and right
            duL = this%grid%elem(i,j)%dudx*rL(1) + this%grid%elem(i,j)%dudy*rL(2)
            duR = this%grid%elem(i,j)%dudx*rR(1) + this%grid%elem(i,j)%dudy*rR(2)

            ! Reconstructing states on left and right of interface
            if (this%limiter=="none") then

              ! Unlimited reconstruction
              do k=1,4
                this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)
                this%grid%edges_v(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + duR(k)
              end do

            else if ((this%limiter=="barth").or.(this%limiter=="venkat")) then

              ! Slope-limited reconstruction
              do k=1,4
                this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + &
                  duL(k)*this%grid%elem(i,j)%phi(k)
                this%grid%edges_v(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + &
                  duR(k)*this%grid%elem(i,j)%phi(k)
              end do

            else
              write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
              stop
            end if

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
            duL = this%grid%elem(i,j)%dudx*rL(1) + this%grid%elem(i,j)%dudy*rL(2)

            ! Reconstructing primitive states on left of interface
            if (this%limiter=="none") then

              ! Unlimited reconstruction
              do k=1,4
                this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)
              end do

            else if ((this%limiter=="barth").or.(this%limiter=="venkat")) then

              ! Slope-limited reconstruction
              do k=1,4
                this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + &
                  duL(k)*this%grid%elem(i,j)%phi(k)
              end do

            else
              write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
              stop
            end if

          else if (j.eq.(this%grid%nelemj)) then

            ! Only need to reconstruct the right state
            rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding the slope on the right side of the interface
            duR = this%grid%elem(i,j)%dudx*rR(1) + this%grid%elem(i,j)%dudy*rR(2)

            ! Reconstructing primitive states on right of interface
            if (this%limiter=="none") then

              ! Unlimited reconstruction
              do k=1,4
                this%grid%edges_h(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + duR(k)
              end do

            else if ((this%limiter=="barth").or.(this%limiter.eq."venkat")) then

              ! Slope-limited reconstruction
              do k=1,4
                this%grid%edges_h(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + &
                  duR(k)*this%grid%elem(i,j)%phi(k)
              end do

            else
              write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
              stop
            end if

          else

            ! Computing position vectors for face midpoints
            rL(1) = this%grid%edges_h(i,j+1)%xm - this%grid%elem(i,j)%xc
            rL(2) = this%grid%edges_h(i,j+1)%ym - this%grid%elem(i,j)%yc
            rR(1) = this%grid%edges_h(i,j)%xm   - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym   - this%grid%elem(i,j)%yc

            ! Finding slopes on left and right
            duL = this%grid%elem(i,j)%dudx*rL(1) + this%grid%elem(i,j)%dudy*rL(2)
            duR = this%grid%elem(i,j)%dudx*rR(1) + this%grid%elem(i,j)%dudy*rR(2)

            ! Reconstructing states on left and right of interface
            if (this%limiter=="none") then
              do k=1,4
                this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)
                this%grid%edges_h(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + duR(k)
              end do
            else if ((this%limiter=="barth").or.(this%limiter=="venkat")) then

              ! Slope-limited reconstruction
              do k=1,4
                this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + &
                  duL(k)*this%grid%elem(i,j)%phi(k)
                this%grid%edges_h(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + &
                  duR(k)*this%grid%elem(i,j)%phi(k)
              end do

            else
              write(*,'(3a)') "Slope limiter ", this%limiter, " not recognized."
              stop
            end if

          end if

        end do
      end do

      ! Must now iterate through all the interior interfaces and solve
      ! the Riemann problem to find the fluxes
      ! Vertical faces
      do j=1,this%grid%nelemj
        do i=2,this%grid%nelemi

          ! Calling Riemann solver
          this%grid%edges_v(i,j)%flux = roe(this%grid%edges_v(i,j)%uL, &
            this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2),this%winfty)

          ! Checking for NaNs
          do k=1,4
            if (isnan(this%grid%edges_v(i,j)%flux(k))) then
              write(*,'(a)') "NaNs encountered after solving Riemann problem for vertical faces"
              write(*,'(a,i4,a,i4,a,i4)') "i = ", i, " j = ", j , " k = ", k
              write(*,'(2(a,f12.5))') "x_midpoint = ", this%grid%edges_v(i,j)%xm, &
                " y_midpoint = ", this%grid%edges_v(i,j)%ym
              stop
            end if
          end do

        end do
      end do

      ! Horizontal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi

          ! Calling Riemann solver
          this%grid%edges_h(i,j)%flux = roe(this%grid%edges_h(i,j)%uL, &
            this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3),this%winfty)

          ! Checking for NaNs
          do k=1,4
            if (isnan(this%grid%edges_h(i,j)%flux(k))) then
              write(*,'(a)') "NaNs encountered after solving Riemann problem for horizontal faces"
              write(*,'(a,i4,a,i4,a,i4)') "i = ", i, " j = ", j , " k = ", k
              write(*,'(2(a,f12.5))') "x_midpoint = ", this%grid%edges_h(i,j)%xm, &
                " y_midpoint = ", this%grid%edges_h(i,j)%ym
              stop
            end if
          end do

        end do
      end do

      ! Applying boundary conditions
      call apply_bcs(this)

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
      double precision                :: duL(4),duR(4)
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
          this%grid%edges_v(i,j)%flux = roe(this%grid%elem(i-1,j)%u, &
            this%grid%elem(i,j)%u,this%grid%elem(i-1,j)%n(:,2),this%winfty)

        end do
      end do

      ! Horizontal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi

          ! Calling Riemann solver
          this%grid%edges_h(i,j)%flux = roe(this%grid%elem(i,j-1)%u, &
            this%grid%elem(i,j)%u,this%grid%elem(i,j-1)%n(:,3),this%winfty)

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
