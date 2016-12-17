!===========================================================
! This module contains a simple solver for the 2D Euler
! equations using a second-order accurate finite volume
! method.
!===========================================================
module solvers
  use cgns
  use mesh_class, only : mesh
  use utils,      only : u_to_w, w_to_u
  use limiters,   only : minmod, vanleer, barth
  use flux,       only : fluxax, fluxay, flux_adv
  use grad,       only : compute_gradient
  use riemann,    only : roe
  implicit none
  private
  public :: solver, initialize, solve_feuler, solve_rk4, residual_inv, write_results_cgns, write_results_tec

  !---------------------------------------------------------
  ! Class for solver
  !---------------------------------------------------------
  type solver
    integer               :: ntsteps            ! Number of time steps
    integer, dimension(4) :: bcids              ! BC identifiers
    double precision      :: dt                 ! Time step
    double precision      :: tfinal             ! Final time
    double precision      :: g                  ! Ratio of specific heats
    double precision      :: winfty(4)          ! Freestream primitive variables
    type(mesh)            :: grid               ! Mesh object
  end type solver

  contains

    !---------------------------------------------------------
    ! Subroutine which initializes solution, i.e., allocates
    ! memory and sets up some important variables
    !---------------------------------------------------------
    subroutine initialize(this,m,delta_t,t_final,gam,w0,winf,bcidents)
      implicit none
      type(solver),       intent(inout)           :: this
      type(mesh),         intent(in)              :: m
      double precision,   intent(in)              :: delta_t,gam,t_final
      double precision,   intent(in), allocatable :: w0(:,:,:)
      double precision,   intent(in)              :: winf(4)
      integer,            intent(in)              :: bcidents(4)
      double precision, allocatable               :: u0(:,:,:)
      integer                                     :: allocate_err,i,j

      print *, "Initializing solver..."

      ! Assigning members of the solver class
      this%grid    = m
      this%dt      = delta_t
      this%tfinal  = t_final
      this%g       = gam
      this%winfty  = winf
      this%ntsteps = ceiling(t_final/delta_t)
      this%bcids   = bcidents
      write (*,'(a,i5)') "number of time steps: ", this%ntsteps

      ! Converting initial condition to conserved variables
      ! (includes ghost cells)
      allocate(u0,mold=w0)
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          u0(i,j,:) = w_to_u(w0(i,j,:),gam)
        end do
      end do

      ! Filling data into elem structures
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%elem(i,j)%u  = u0(i,j,:)
          this%grid%elem(i,j)%w  = w0(i,j,:)
        end do
      end do

      print *, "Done initializing solver."

    end subroutine initialize

    !---------------------------------------------------------
    ! Subroutine for solving to some final time using
    ! forward Euler time stepping
    !---------------------------------------------------------
    subroutine solve_feuler(this,write_freq)
      implicit none
      type(solver), intent(inout) :: this
      integer, intent(in) :: write_freq
      double precision, allocatable :: u(:,:,:), resid(:,:,:)
      character (len=30) :: tecname
      integer :: i,j,k,aer

      ! Allocating memory for the solution and the residual
      allocate(resid(this%grid%nelemi,this%grid%nelemj,4),stat=aer)
      if (aer.ne.0) then
        print *, "Error: can't allocate memory for resid in solve_feuler."
        stop
      end if

      ! Writing initial solution
      write (tecname, '(a,i0,a)') "sol", 0, ".tec"
      print *, "Writing data to ", tecname
      call write_results_tec(this,tecname)

      ! Marching in time
      do k=1,this%ntsteps

        ! Printing some information
        write(*,'(a,i5,a,es12.5)') "timestep: ", k, " t = ", dble(k)*this%dt

        ! Writing current solution
        if (mod(k,write_freq).eq.0) then
          write (tecname, '(a,i0,a)') "sol", (k), ".tec"
          print *, "Writing data to ", tecname
          call write_results_tec(this,tecname)
        end if

        ! Copying old solution
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%u0 = this%grid%elem(i,j)%u
          end do
        end do

        ! Computing residual
        call residual_inv(this,resid)

        ! Advancing in time
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt*resid(i,j,:)
          end do
        end do

      end do

    end subroutine solve_feuler


    !---------------------------------------------------------
    ! Subroutine for solving to some final time using
    ! Runge-Kutta 4
    !---------------------------------------------------------
    subroutine solve_rk4(this,write_freq)
      implicit none
      type(solver), intent(inout) :: this
      integer, intent(in) :: write_freq
      double precision, allocatable :: u(:,:,:), resid(:,:,:)
      character (len=30) :: tecname
      integer :: i,j,k,aer

      ! Allocating memory for the solution and the residual
      allocate(resid(this%grid%nelemi,this%grid%nelemj,4),stat=aer)
      if (aer.ne.0) then
        print *, "Error: can't allocate memory for resid in solve_feuler."
        stop
      end if

      ! Writing initial solution
      write (tecname, '(a,i0,a)') "sol", 0, ".tec"
      print *, "Writing data to ", tecname
      call write_results_tec(this,tecname)

      ! Marching in time
      do k=1,this%ntsteps

        ! Printing some information
        write(*,'(a,i5,a,es12.5)') "timestep: ", k, " t = ", dble(k)*this%dt

        ! Writing current solution
        if (mod(k,write_freq).eq.0) then
          write (tecname, '(a,i0,a)') "sol", (k), ".tec"
          print *, "Writing data to ", tecname
          call write_results_tec(this,tecname)
        end if

        ! Copying old solution
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%u0 = this%grid%elem(i,j)%u
          end do
        end do

        ! Stage 1
        call residual_inv(this,resid)
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt/4.0*resid(i,j,:)
          end do
        end do

        ! Stage 2
        call residual_inv(this,resid)
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt/3.0*resid(i,j,:)
          end do
        end do

        ! Stage 3
        call residual_inv(this,resid)
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt/2.0*resid(i,j,:)
          end do
        end do

        ! Stage 4
        call residual_inv(this,resid)
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt*resid(i,j,:)
          end do
        end do

      end do

    end subroutine solve_rk4

    !---------------------------------------------------------
    ! Subroutine for applying boundary conditions
    ! 1000 - farfield     (primitives specified)
    ! 1001 - extrapolate  (no information needed)
    ! 1002 - slip wall    (no information needed)
    ! 1003 - no-slip wall (wall temperature, assumed adiabatic)
    !---------------------------------------------------------
    subroutine apply_bcs(this)
      implicit none
      type(solver), intent(inout)    :: this
      double precision, dimension(4) :: wextrap,uextrap,wtmp
      double precision               :: pw
      !double precision               :: m,p1,p2,d1,d2
      integer                        :: i,j

      !===============================================
      ! Bottom boundary (inward pointing normal)
      !===============================================
      j = 1
      if (this%bcids(1).eq.1000) then

        ! Farfield
        do i=1,this%grid%nelemi
          this%grid%edges_h(i,j)%flux = -flux_adv(this%winfty,this%grid%elem(i,j)%n(:,1),this%g)
        end do

      else if (this%bcids(1).eq.1001) then

        ! Extrapolate
        do i=1,this%grid%nelemi
          uextrap = this%grid%elem(i,j)%u
          wextrap = u_to_w(uextrap,this%g)
          this%grid%edges_h(i,j)%flux = -flux_adv(wextrap,this%grid%elem(i,j)%n(:,1),this%g)
        end do

      else if (this%bcids(1).eq.1002) then

        ! Slip wall
        do i=1,this%grid%nelemi
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = wtmp(4)
          this%grid%edges_h(i,j)%flux(1) = 0.0d0
          this%grid%edges_h(i,j)%flux(2) = -pw*this%grid%elem(i,j)%n(1,1)
          this%grid%edges_h(i,j)%flux(3) = -pw*this%grid%elem(i,j)%n(2,1)
          this%grid%edges_h(i,j)%flux(4) = 0.0d0
        end do

      else if (this%bcids(1).eq.1003) then

        ! No-slip wall
        print *, "Warning: no-slip wall bc not implemented yet."

      else

        print *, "Boundary condition ", this%bcids(1), " for bottom not recognized."
        stop

      end if

      !===============================================
      ! Right boundary
      !===============================================
      i = this%grid%nelemi
      if (this%bcids(2).eq.1000) then

        ! Farfield
        do j=1,this%grid%nelemj
          this%grid%edges_v(i+1,j)%flux = flux_adv(this%winfty,this%grid%elem(i,j)%n(:,2),this%g)
        end do

      else if (this%bcids(2).eq.1001) then

        ! Extrapolate
        do j=1,this%grid%nelemj
          uextrap = this%grid%elem(i,j)%u
          wextrap = u_to_w(uextrap,this%g)
          this%grid%edges_v(i+1,j)%flux = flux_adv(wextrap,this%grid%elem(i,j)%n(:,2),this%g)
        end do

      else if (this%bcids(2).eq.1002) then

        ! Slip wall
        do i=1,this%grid%nelemi
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = wtmp(4)
          this%grid%edges_v(i+1,j)%flux(1) = 0.0d0
          this%grid%edges_v(i+1,j)%flux(2) = pw*this%grid%elem(i,j)%n(1,2)
          this%grid%edges_v(i+1,j)%flux(3) = pw*this%grid%elem(i,j)%n(2,2)
          this%grid%edges_v(i+1,j)%flux(4) = 0.0d0
        end do

      else if (this%bcids(2).eq.1003) then

        ! No-slip wall
        print *, "Warning: no-slip wall bc not implemented yet."

      else

        print *, "Boundary condition ", this%bcids(2), " for right not recognized."
        stop

      end if

      !===============================================
      ! Top boundary
      !===============================================
      j = this%grid%nelemj
      if (this%bcids(3).eq.1000) then

        ! Farfield
        do i=1,this%grid%nelemi
          this%grid%edges_h(i,j+1)%flux = flux_adv(this%winfty,this%grid%elem(i,j)%n(:,3),this%g)
        end do

      else if (this%bcids(3).eq.1001) then

        ! Extrapolate
        do i=1,this%grid%nelemi
          uextrap = this%grid%elem(i,j)%u
          wextrap = u_to_w(uextrap,this%g)
          this%grid%edges_h(i,j+1)%flux = flux_adv(wextrap,this%grid%elem(i,j)%n(:,3),this%g)
        end do

      else if (this%bcids(3).eq.1002) then

        ! Slip wall
        do i=1,this%grid%nelemi
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = wtmp(4)
          this%grid%edges_h(i,j+1)%flux(1) = 0.0d0
          this%grid%edges_h(i,j+1)%flux(2) = pw*this%grid%elem(i,j)%n(1,3)
          this%grid%edges_h(i,j+1)%flux(3) = pw*this%grid%elem(i,j)%n(2,3)
          this%grid%edges_h(i,j+1)%flux(4) = 0.0d0
        end do

      else if (this%bcids(3).eq.1003) then

        ! No-slip wall
        print *, "Warning: no-slip wall bc not implemented yet."

      else

        print *, "Boundary condition ", this%bcids(3), " for top not recognized."
        stop

      end if

      !===============================================
      ! Left boundary (inward pointing normal)
      !===============================================
      i = 1
      if (this%bcids(4).eq.1000) then

        ! Farfield
        do j=1,this%grid%nelemj
          this%grid%edges_v(i,j)%flux = -flux_adv(this%winfty,this%grid%elem(i,j)%n(:,4),this%g)
        end do

      else if (this%bcids(4).eq.1001) then

        ! Extrapolate
        do j=1,this%grid%nelemj
          uextrap = this%grid%elem(i,j)%u
          wextrap = u_to_w(uextrap,this%g)
          this%grid%edges_v(i,j)%flux = -flux_adv(wextrap,this%grid%elem(i,j)%n(:,4),this%g)
        end do

      else if (this%bcids(4).eq.1002) then

        ! Slip wall
        do i=1,this%grid%nelemi
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = wtmp(4)
          this%grid%edges_v(i,j)%flux(1) = 0.0d0
          this%grid%edges_v(i,j)%flux(2) = -pw*this%grid%elem(i,j)%n(1,4)
          this%grid%edges_v(i,j)%flux(3) = -pw*this%grid%elem(i,j)%n(2,4)
          this%grid%edges_v(i,j)%flux(4) = 0.0d0
        end do

      else if (this%bcids(4).eq.1003) then

        ! No-slip wall
        print *, "Warning: no-slip wall bc not implemented yet."

      else

        print *, "Boundary condition ", this%bcids(4), " for left not recognized."
        stop

      end if

    end subroutine apply_bcs

    !---------------------------------------------------------
    ! Subroutine for computing the residual
    !---------------------------------------------------------
    subroutine residual_inv(this,resid)
      implicit none
      type(solver), intent(inout)     :: this
      double precision, intent(inout) :: resid(:,:,:)
      double precision                :: duL(4),duR(4)
      double precision, dimension(4)  :: fx,fy,uextrap,wextrap,wtmp
      double precision, dimension(2)  :: rL,rR,r
      double precision, allocatable   :: u(:,:,:),gradU(:,:,:)
      integer                         :: i,j,k,err

      ! Allocating memory for gradient of state at cell centers
      allocate(gradU(this%grid%nelemi,this%grid%nelemj,8),stat=err)
      if (err.ne.0) then
        print *, "Error: Can't allocate memory for gradU in update."
        stop
      end if
      allocate(u(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for u in solve_feuler."
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

      ! Reconstructing states on left and right sides of each vertical interface
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi

          if (i.eq.1) then

            ! Only need to reconstruct the left state
            rL(1) = this%grid%edges_v(i+1,j)%xm - this%grid%elem(i,j)%xc
            rL(2) = this%grid%edges_v(i+1,j)%ym - this%grid%elem(i,j)%yc

            ! Finding the slope on the left side of the interface
            duL = gradU(i,j,1:4)*rL(1) + gradU(i,j,5:8)*rL(2)

            ! Reconstructing state on left of interface
            do k=1,4

              ! Reconstruction
              !this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + barth(duL(k),this%grid%elem(i,j)%u(k),umax,umin)
              this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)

            end do

          else if (i.eq.(this%grid%nelemi)) then

            ! Only need to reconstruct the right state
            rR(1) = this%grid%edges_v(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_v(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding the slope on the right side of the interface
            duR = gradU(i,j,1:4)*rR(1) + gradU(i,j,5:8)*rR(2)

            ! Reconstructing primitive states on right of interface
            do k=1,4

              ! Reconstruction
              !this%grid%edges_v(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + barth(duR(k),this%grid%elem(i,j)%u(k),umax,umin)
              this%grid%edges_v(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + duR(k)

            end do
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
            do k=1,4

              ! Unlimited
              this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)
              this%grid%edges_v(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + duR(k)

              ! Not sure that this is right so I'm going to go with the Barth limiter
              !this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)*minmod(r(k))
              !this%grid%edges_v(i+1,j)%uR(k) = this%grid%elem(i,j)%u(k) + duR(k)*minmod(1.0d0/r(k))

            end do

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
            do k=1,4

              ! Reconstruction
              !this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + barth(duL(k),this%grid%elem(i,j)%u(k),umax,umin)
              this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)

            end do

          else if (j.eq.(this%grid%nelemj)) then

            ! Only need to reconstruct the right state
            rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding the slope on the right side of the interface
            duR = gradU(i,j,1:4)*rR(1) + gradU(i,j,5:8)*rR(2)

            ! Reconstructing primitive states on right of interface
            do k=1,4

              ! Reconstruction
              !this%grid%edges_h(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + barth(duR(k),this%grid%elem(i,j)%u(k),umax,umin)
              this%grid%edges_h(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + duR(k)

            end do
          else

            ! Computing position vectors for face midpoints
            rL(1) = this%grid%edges_h(i,j+1)%xm - this%grid%elem(i,j)%xc
            rL(2) = this%grid%edges_h(i,j+1)%ym - this%grid%elem(i,j)%yc
            rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding slopes on left and right
            duL = gradU(i,j,1:4)*rL(1) + gradU(i,j,5:8)*rL(2)
            duR = gradU(i,j,1:4)*rR(1) + gradU(i,j,5:8)*rR(2)

            ! Reconstructing primitive states on left and right of interface
            do k=1,4

              this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + duL(k)
              this%grid%edges_h(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + duR(k)

            end do

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
          this%grid%edges_v(i,j)%flux = roe(this%grid%edges_v(i,j)%uL,this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2))

        end do
      end do
      ! Horizontal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi

          ! Calling Riemann solver
          this%grid%edges_h(i,j)%flux = roe(this%grid%edges_h(i,j)%uL,this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3))

        end do
      end do

      ! Advance one step in time (forward Euler for now)
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

    !---------------------------------------------------------
    ! Subroutine for writing results out to a CGNS file
    !---------------------------------------------------------
    subroutine write_results_cgns(this, file_name)
      implicit none
      type(solver), intent(inout)  :: this
      character (len=*),intent(in) :: file_name
      integer (kind=4)             :: isize(2,3)
      character (len=30)           :: basename,zonename
      integer                      :: i,j,ier
      integer                      :: idx_file,idx_base
      integer                      :: icelldim,iphysdim
      integer, dimension(3)        :: irmin,irmax

      ! Setting inputs for CGNS file
      idx_file   = 1
      idx_base   = 1
      icelldim   = 2
      iphysdim   = 2
      isize(1,1) = this%grid%imax
      isize(2,1) = this%grid%jmax
      isize(1,2) = this%grid%nelemi
      isize(2,2) = this%grid%nelemj
      isize(1,3) = 0
      isize(2,3) = 0
      basename   = "Base"
      zonename   = "Solution"

      ! Opening CGNS file
      call cg_open_f(file_name,CG_MODE_WRITE,idx_file,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

      ! Creating zone

      ! Closing CGNS file
      call cg_close_f(idx_file,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

    end subroutine write_results_cgns

    !---------------------------------------------------------
    ! Subroutine for writing results out to a Tecplot file
    !---------------------------------------------------------
    subroutine write_results_tec(this, file_name)
      implicit none
      type(solver), intent(in)  :: this
      character (len=*),intent(in) :: file_name
      integer ::i,j,k

        ! Writing result to file
        open(2,file=file_name)
        write(2,'(a,i5,a)') 'title="Step ', (k-1), '"'
        write(2,'(a)') 'variables="x","y","rho","rhou","rhov","E"'
        write(2,'(a,i5,a,i5)') 'zone i=', this%grid%imax, ' j=', this%grid%jmax
        write(2,'(a)') 'datapacking=block'
        write(2,'(a)') 'varlocation=([3,4,5,6]=cellcentered)'
        do j=1,this%grid%jmax
          do i=1,this%grid%imax
            write(2,'(es25.10)',advance='no') this%grid%x(i,j)
          end do
          write(2,'(a)') " "
        end do
        do j=1,this%grid%jmax
          do i=1,this%grid%imax
            write(2,'(es25.10)',advance='no') this%grid%y(i,j)
          end do
          write(2,'(a)') " "
        end do
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            write(2,'(es25.10)',advance='no') this%grid%elem(i,j)%u(1)
          end do
          write(2,'(a)') " "
        end do
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            write(2,'(es25.10)',advance='no') this%grid%elem(i,j)%u(2)
          end do
          write(2,'(a)') " "
        end do
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            write(2,'(es25.10)',advance='no') this%grid%elem(i,j)%u(3)
          end do
          write(2,'(a)') " "
        end do
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            write(2,'(es25.10)',advance='no') this%grid%elem(i,j)%u(4)
          end do
          write(2,'(a)') " "
        end do
        close(2)

    end subroutine write_results_tec

end module solvers
