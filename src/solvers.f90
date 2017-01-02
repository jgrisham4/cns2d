!===========================================================
! This module contains a simple solver for the 2D Euler
! equations and Navier-Stokes equations on curvilinear
! meshes using an upwind second-order accurate finite
! volume method. The discretization makes use of Roe flux
! difference splitting, piecewise linear reconstruction,
! the barth slope limiter, first-order accuracy in time
! using the forward Euler method or fourth-order accuracy
! in time using a Runge Kutta 4 method.
!
! Author: James Grisham
! Date: 01-01-2017
!===========================================================

module solvers
  use cgns
  use mesh_class,      only : mesh
  use utils,           only : u_to_w, w_to_u, compute_elem_max, compute_elem_min
  use limiters,        only : minmod, vanleer
  use flux,            only : fluxax, fluxay, flux_adv
  use grad,            only : compute_gradient
  use riemann,         only : roe,rotated_rhll
  use mms,             only : s_continuity,s_xmom,s_ymom,s_energy,rho_e,u_e,v_e,et_e
  use ieee_arithmetic, only : ieee_is_finite
  implicit none
  private :: apply_bcs
  public  :: solver, initialize, solve_feuler, solve_rk4, residual_inv, write_results_cgns, write_results_tec

  !---------------------------------------------------------
  ! Class for solver
  !---------------------------------------------------------
  type solver
    integer               :: ntsteps   ! Number of time steps
    integer, dimension(4) :: bcids     ! BC identifiers
    double precision      :: dt        ! Time step
    double precision      :: tfinal    ! Final time
    double precision      :: g         ! Ratio of specific heats
    double precision      :: R         ! Ideal gas constant for air
    double precision      :: winfty(4) ! Freestream primitive variables
    character (len=30)    :: limiter   ! Name of slope limiter ("none" or "barth")
    type(mesh)            :: grid      ! Mesh object
  end type solver

  contains

    !---------------------------------------------------------
    ! Subroutine which initializes solution, i.e., allocates
    ! memory and sets up some important variables
    !---------------------------------------------------------
    subroutine initialize(this,m,delta_t,t_final,gam,R,w0,winf,bcidents,lim)
      implicit none
      type(solver),       intent(inout)           :: this
      type(mesh),         intent(in)              :: m
      double precision,   intent(in)              :: delta_t,gam,t_final,R
      double precision,   intent(in), allocatable :: w0(:,:,:)
      double precision,   intent(in)              :: winf(4)
      integer,            intent(in)              :: bcidents(4)
      character (len=*),  intent(in)              :: lim
      double precision, allocatable               :: u0(:,:,:)
      integer                                     :: allocate_err,i,j

      print *, "Initializing solver..."

      ! Assigning members of the solver class
      this%grid    = m
      this%dt      = delta_t
      this%tfinal  = t_final
      this%g       = gam
      this%R       = R
      this%winfty  = winf
      this%ntsteps = nint(t_final/delta_t)
      this%bcids   = bcidents
      this%limiter = lim
      write (*,'(a,i7)') "number of time steps: ", this%ntsteps

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
            this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt/4.0d0*resid(i,j,:)
          end do
        end do

        ! Stage 2
        call residual_inv(this,resid)
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt/3.0d0*resid(i,j,:)
          end do
        end do

        ! Stage 3
        call residual_inv(this,resid)
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt/2.0d0*resid(i,j,:)
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
    ! 1000 - farfield         (primitives specified)
    ! 1001 - extrapolate      (no information needed)
    ! 1002 - slip wall weak   (no information needed)
    ! 1003 - slip wall strong (no information needed)
    ! 1004 - no-slip wall     (wall temperature, assumed adiabatic)
    ! 2000 - mms solution     (Dirichlet along walls--state is set)
    !---------------------------------------------------------
    subroutine apply_bcs(this)
      implicit none
      type(solver), intent(inout)    :: this
      double precision, dimension(4) :: wextrap,uextrap,utmp,wtmp
      double precision, dimension(4) :: duds
      double precision               :: pw,ds,rL(2),rR(2)
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

        ! Weak enforcement of slip wall
        ! This is accomplished by setting the fluxes using
        ! extrapolated pressure at the wall.
        do i=1,this%grid%nelemi
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = 3.0d0/2.0d0*wtmp(4)
          wtmp = u_to_w(this%grid%elem(i,j+1)%u,this%g)
          pw = pw - 0.5d0*wtmp(4)
          this%grid%edges_h(i,j)%flux(1) = 0.0d0
          this%grid%edges_h(i,j)%flux(2) = -pw*this%grid%elem(i,j)%n(1,1)
          this%grid%edges_h(i,j)%flux(3) = -pw*this%grid%elem(i,j)%n(2,1)
          this%grid%edges_h(i,j)%flux(4) = 0.0d0
        end do

      else if (this%bcids(1).eq.1003) then

        ! Strong enforcement of slip wall
        ! This is accomplished by setting the state at the ghost cells
        ! Reconstructing the interface states and solving the Riemann
        ! problem at the interface
        do i=1,this%grid%nelemi

          ! Setting state in the ghost cells
          this%grid%elem(i,j-1)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i,j-1)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i,j-1)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i,j-1)%u(4) =  this%grid%elem(i,j)%u(4)
          this%grid%elem(i,j-2)%u(1) =  this%grid%elem(i,j+1)%u(1)
          this%grid%elem(i,j-2)%u(2) = -this%grid%elem(i,j+1)%u(2)
          this%grid%elem(i,j-2)%u(3) = -this%grid%elem(i,j+1)%u(3)
          this%grid%elem(i,j-2)%u(4) =  this%grid%elem(i,j+1)%u(4)

          ! Finding the gradient using first-order accurate differences
          ds = sqrt((this%grid%elem(i,j-1)%xc - this%grid%elem(i,j-2)%xc)**2 + &
                    (this%grid%elem(i,j-1)%yc - this%grid%elem(i,j-2)%yc)**2)
          duds = (this%grid%elem(i,j-1)%u - this%grid%elem(i,j-2)%u)/ds

          ! Computing the left state of the interface
          rL(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j-1)%xc
          rL(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j-1)%yc
          this%grid%edges_h(i,j)%uL = this%grid%elem(i,j-1)%u + duds*sqrt(rL(1)**2 + rL(2)**2)

          ! Computing the right state of the interface
          rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
          rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc
          this%grid%edges_h(i,j)%uR = this%grid%elem(i,j)%u + &
                                      this%grid%elem(i,j)%dudx*rR(1) + &
                                      this%grid%elem(i,j)%dudy*rR(2)

          ! Solving the Riemann problem at the interface
          this%grid%edges_h(i,j)%flux = roe(this%grid%edges_h(i,j)%uL,this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3))

        end do

      else if (this%bcids(1).eq.2000) then

        do i=1,this%grid%nelemi

          ! Finding the conserved variables from the exact solution
          utmp(1) = rho_e(this%grid%edges_h(i,j)%xm,this%grid%edges_h(i,j)%ym)
          utmp(2) = utmp(1)*u_e(this%grid%edges_h(i,j)%xm,this%grid%edges_h(i,j)%ym)
          utmp(3) = utmp(1)*v_e(this%grid%edges_h(i,j)%xm,this%grid%edges_h(i,j)%ym)
          utmp(4) = et_e(this%grid%edges_h(i,j)%xm,this%grid%edges_h(i,j)%ym)

          ! I MUST ADD ALL THE VISCOUS FLUXES TO THE BOUNDARY CONDITIONS

          ! Finding the primitive variables at the face and computing the flux
          wtmp = u_to_w(utmp,this%g)
          this%grid%edges_h(i,j)%flux = -(flux_adv(wtmp,this%grid%elem(i,j)%n(:,1),this%g) + &
            flux_visc(wtmp,)

        end do

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
        do j=1,this%grid%nelemj
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = 3.0d0/2.0d0*wtmp(4)
          wtmp = u_to_w(this%grid%elem(i-1,j)%u,this%g)
          pw = pw - 0.5d0*wtmp(4)
          this%grid%edges_v(i+1,j)%flux(1) = 0.0d0
          this%grid%edges_v(i+1,j)%flux(2) = pw*this%grid%elem(i,j)%n(1,2)
          this%grid%edges_v(i+1,j)%flux(3) = pw*this%grid%elem(i,j)%n(2,2)
          this%grid%edges_v(i+1,j)%flux(4) = 0.0d0
        end do

      else if (this%bcids(2).eq.1003) then

        ! Slip wall strongly enforced
        do j=1,this%grid%nelemj

          ! Setting state in ghost cells
          this%grid%elem(i+1,j)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i+1,j)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i+1,j)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i+1,j)%u(4) =  this%grid%elem(i,j)%u(4)
          this%grid%elem(i+2,j)%u(1) =  this%grid%elem(i-1,j)%u(1)
          this%grid%elem(i+2,j)%u(2) = -this%grid%elem(i-1,j)%u(2)
          this%grid%elem(i+2,j)%u(3) = -this%grid%elem(i-1,j)%u(3)
          this%grid%elem(i+2,j)%u(4) =  this%grid%elem(i-1,j)%u(4)

          ! Computing the gradient in the ghost cells using
          ! first-order accurate differences
          ds = sqrt((this%grid%edges_v(i+1,j)%xm - this%grid%elem(i+1,j)%xc)**2 + &
                    (this%grid%edges_v(i+1,j)%ym - this%grid%elem(i+1,j)%yc)**2)
          duds = (this%grid%elem(i+1,j)%u - this%grid%elem(i+2,j)%u)/ds

          ! Reconstructing the state on the interior
          rL(1) = this%grid%edges_v(i+1,j)%xm - this%grid%elem(i,j)%xc
          rL(2) = this%grid%edges_v(i+1,j)%ym - this%grid%elem(i,j)%yc
          this%grid%edges_v(i+1,j)%uL = this%grid%elem(i,j)%u + this%grid%elem(i,j)%dudx*rL(1) + &
            this%grid%elem(i,j)%dudy*rL(2)

          ! Reconstructing the state on the exterior
          rR(1) = this%grid%edges_v(i+1,j)%xm - this%grid%elem(i+1,j)%xc
          rR(2) = this%grid%edges_v(i+1,j)%ym - this%grid%elem(i+1,j)%yc
          this%grid%edges_v(i+1,j)%uR = this%grid%elem(i+1,j)%u + duds*sqrt(rR(1)**2+rR(2)**2)

          ! Solving the Riemann problem at the interface
          this%grid%edges_v(i+1,j)%flux = roe(this%grid%edges_v(i+1,j)%uL,this%grid%edges_v(i+1,j)%uR,this%grid%elem(i,j)%n(:,2))

        end do

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
          pw = 3.0d0/2.0d0*wtmp(4)
          wtmp = u_to_w(this%grid%elem(i,j-1)%u,this%g)
          pw = pw - 0.5d0*wtmp(4)
          this%grid%edges_h(i,j+1)%flux(1) = 0.0d0
          this%grid%edges_h(i,j+1)%flux(2) = pw*this%grid%elem(i,j)%n(1,3)
          this%grid%edges_h(i,j+1)%flux(3) = pw*this%grid%elem(i,j)%n(2,3)
          this%grid%edges_h(i,j+1)%flux(4) = 0.0d0
        end do

      else if (this%bcids(3).eq.1003) then

        ! Strong enforcement of slip wall
        ! This is accomplished by setting the state at the ghost cells
        ! Reconstructing the interface states and solving the Riemann
        ! problem at the interface
        do i=1,this%grid%nelemi

          ! Setting state in the ghost cells
          this%grid%elem(i,j+1)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i,j+1)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i,j+1)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i,j+1)%u(4) =  this%grid%elem(i,j)%u(4)
          this%grid%elem(i,j+2)%u(1) =  this%grid%elem(i,j-1)%u(1)
          this%grid%elem(i,j+2)%u(2) = -this%grid%elem(i,j-1)%u(2)
          this%grid%elem(i,j+2)%u(3) = -this%grid%elem(i,j-1)%u(3)
          this%grid%elem(i,j+2)%u(4) =  this%grid%elem(i,j-1)%u(4)

          ! Finding the gradient using first-order accurate differences
          ds = sqrt((this%grid%elem(i,j+1)%xc - this%grid%elem(i,j+2)%xc)**2 + &
                    (this%grid%elem(i,j+1)%yc - this%grid%elem(i,j+2)%yc)**2)
          duds = (this%grid%elem(i,j+1)%u - this%grid%elem(i,j+2)%u)/ds

          ! o------------o
          ! |            | j+2 -- ghost
          ! |            |
          ! o------------o
          ! |            | j+1 -- ghost
          ! | R          |
          ! o------------o
          ! | L          | j
          ! |            |
          ! o------------o
          ! Computing the left state of the interface
          rL(1) = this%grid%edges_h(i,j+1)%xm - this%grid%elem(i,j)%xc
          rL(2) = this%grid%edges_h(i,j+1)%ym - this%grid%elem(i,j)%yc
          this%grid%edges_h(i,j+1)%uL = this%grid%elem(i,j)%u + &
                                        this%grid%elem(i,j)%dudx*rL(1) + &
                                        this%grid%elem(i,j)%dudy*rL(2)

          ! Computing the right state of the interface
          rR(1) = this%grid%edges_h(i,j+1)%xm - this%grid%elem(i,j+1)%xc
          rR(2) = this%grid%edges_h(i,j+1)%ym - this%grid%elem(i,j+1)%yc
          this%grid%edges_h(i,j+1)%uR = this%grid%elem(i,j+1)%u + duds*sqrt(rR(1)**2 + rR(2)**2)

          ! Solving the Riemann problem at the interface
          this%grid%edges_h(i,j+1)%flux = roe(this%grid%edges_h(i,j+1)%uL,this%grid%edges_h(i,j+1)%uR,this%grid%elem(i,j)%n(:,3))

        end do

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
        do j=1,this%grid%nelemj
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = 3.0d0/2.0d0*wtmp(4)
          wtmp = u_to_w(this%grid%elem(i+1,j)%u,this%g)
          pw = pw - 0.5d0*wtmp(4)
          this%grid%edges_v(i,j)%flux(1) = 0.0d0
          this%grid%edges_v(i,j)%flux(2) = -pw*this%grid%elem(i,j)%n(1,4)
          this%grid%edges_v(i,j)%flux(3) = -pw*this%grid%elem(i,j)%n(2,4)
          this%grid%edges_v(i,j)%flux(4) = 0.0d0
        end do

      else if (this%bcids(4).eq.1003) then

        ! Strongly enforced slip wall
        do j=1,this%grid%nelemj

          ! Assigning states in the ghost cells
          this%grid%elem(i-1,j)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i-1,j)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i-1,j)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i-1,j)%u(4) =  this%grid%elem(i,j)%u(4)
          this%grid%elem(i-2,j)%u(1) =  this%grid%elem(i+1,j)%u(1)
          this%grid%elem(i-2,j)%u(2) = -this%grid%elem(i+1,j)%u(2)
          this%grid%elem(i-2,j)%u(3) = -this%grid%elem(i+1,j)%u(3)
          this%grid%elem(i-2,j)%u(4) =  this%grid%elem(i+1,j)%u(4)

          ! Computing gradient in the ghost cell
          ds = sqrt((this%grid%elem(i-1,j)%xc - this%grid%elem(i-2,j)%xc)**2 + &
                    (this%grid%elem(i-1,j)%yc - this%grid%elem(i-2,j)%yc)**2)
          duds = (this%grid%elem(i-1,j)%u - this%grid%elem(i-2,j)%u)/ds

          ! Reconstructing the exterior state
          rL(1) = this%grid%edges_v(i,j)%xm - this%grid%elem(i-1,j)%xc
          rL(2) = this%grid%edges_v(i,j)%ym - this%grid%elem(i-1,j)%yc
          this%grid%edges_v(i,j)%uL = this%grid%elem(i-1,j)%u + duds*sqrt(rL(1)**2 + rL(2)**2)

          ! Reconstructing the interior state
          rR(1) = this%grid%edges_v(i,j)%xm - this%grid%elem(i,j)%xc
          rR(2) = this%grid%edges_v(i,j)%ym - this%grid%elem(i,j)%yc
          this%grid%edges_v(i,j)%uR = this%grid%elem(i,j)%u + &
                                      this%grid%elem(i,j)%dudx*rR(1) + &
                                      this%grid%elem(i,j)%dudy*rR(2)

          ! Solving the Riemann problem
          this%grid%edges_v(i,j)%flux = roe(this%grid%edges_v(i,j)%uL,this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2))

        end do

      else

        print *, "Boundary condition ", this%bcids(4), " for left not recognized."
        stop

      end if

    end subroutine apply_bcs

    !---------------------------------------------------------
    ! Subroutine for computing the residual for Euler eqs
    !---------------------------------------------------------
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
            rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc

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
          this%grid%edges_v(i,j)%flux = roe(this%grid%edges_v(i,j)%uL,this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2))
          !this%grid%edges_v(i,j)%flux = rotated_rhll(this%grid%edges_v(i,j)%uL,this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2))

          ! Checking for NaNs
          do k=1,4
            if (isnan(this%grid%edges_v(i,j)%flux(k))) then
              write(*,'(a)') "NaNs encountered after solving Riemann problem for vertical faces"
              write(*,'(a,i4,a,i4,a,i4)') "i = ", i, " j = ", j , " k = ", k
              stop
            end if
          end do

        end do
      end do

      ! Horizontal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi

          ! Calling Riemann solver
          this%grid%edges_h(i,j)%flux = roe(this%grid%edges_h(i,j)%uL,this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3))
          !this%grid%edges_h(i,j)%flux = rotated_rhll(this%grid%edges_h(i,j)%uL,this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3))

          ! Checking for NaNs
          do k=1,4
            if (isnan(this%grid%edges_h(i,j)%flux(k))) then
              write(*,'(a)') "NaNs encountered after solving Riemann problem for horizontal faces"
              write(*,'(a,i4,a,i4,a,i4)') "i = ", i, " j = ", j , " k = ", k
              stop
            end if
          end do

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

    !---------------------------------------------------------
    ! Subroutine for computing the residual for Navier-Stokes
    ! equations
    !---------------------------------------------------------
    subroutine residual_visc(this,resid)
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
        print *, "Error: Can't allocate memory for gradU in residual_visc."
        stop
      end if
      allocate(u(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for u in residual_visc."
        stop
      end if
      allocate(umax(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for umax in residual_visc."
        stop
      end if
      allocate(umin(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for umin in residual_visc."
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
            rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc

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
          this%grid%edges_v(i,j)%flux = roe(this%grid%edges_v(i,j)%uL,this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2))
          !this%grid%edges_v(i,j)%flux = rotated_rhll(this%grid%edges_v(i,j)%uL,this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2))

          ! Checking for NaNs
          do k=1,4
            if (isnan(this%grid%edges_v(i,j)%flux(k))) then
              write(*,'(a)') "NaNs encountered after solving Riemann problem for vertical faces"
              write(*,'(a,i4,a,i4,a,i4)') "i = ", i, " j = ", j , " k = ", k
              stop
            end if
          end do

        end do
      end do

      ! Horizontal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi

          ! Calling Riemann solver
          this%grid%edges_h(i,j)%flux = roe(this%grid%edges_h(i,j)%uL,this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3))
          !this%grid%edges_h(i,j)%flux = rotated_rhll(this%grid%edges_h(i,j)%uL,this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3))

          ! Checking for NaNs
          do k=1,4
            if (isnan(this%grid%edges_h(i,j)%flux(k))) then
              write(*,'(a)') "NaNs encountered after solving Riemann problem for horizontal faces"
              write(*,'(a,i4,a,i4,a,i4)') "i = ", i, " j = ", j , " k = ", k
              stop
            end if
          end do

        end do
      end do

      ! Computing gradient of temperature
      ! dT/dx_ij = gradU(i,j,1)
      ! dT/dy_ij = gradU(i,j,5)
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%elem(i,j)%w = u_to_w(this%grid%elem(i,j)%u,this%g)
          u(i,j,1) = this%grid%elem(i,j)%w(4)/(this%R*this%grid%elem(i,j)%w(1))
          u(i,j,2:4) = 0.0d0
        end do
      end do
      call compute_gradient(this%grid,u,gradU)

      ! Copying gradient of temperature to element objects
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%elem(i,j)%dTdx = gradU(i,j,1)
          this%grid%elem(i,j)%dTdy = gradU(i,j,5)
        end do
      end do

      ! Computing viscous fluxes for vertical internal faces
      do j=1,this%grid%nelemj
        do i=2,this%grid%nelemi
          this%grid%edges_v(i,j)%flux = this%grid%edges_v(i,j)%flux - &
            flux_visc(this%grid%elem(i-1,j),this%grid%elem(i,j), &
                      this%grid%edges_v(i,j)%n(:,2),this%g,this%R)
        end do
      end do

      ! Computing viscous fluxes for the horizontal internal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%edges_h(i,j)%flux = this%grid%edges_h(i,j)%flux - &
            flux_visc(this%grid%elem(i-1,j),this%grid%elem(i,j), &
                      this%grid%edges_v(i,j)%n(:,3),this%g,this%R)
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

    end subroutine residual_visc

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

    !-------------------------------------------------------
    ! Subroutine for the Barth-Jespersen limiter
    !-------------------------------------------------------
    pure function barth(s,i,j,gradu,umax,umin) result(ph)
      implicit none
      type(solver), intent(in)                  :: s
      integer, intent(in)                       :: i,j
      double precision, allocatable, intent(in) :: gradu(:,:,:),umax(:,:,:),umin(:,:,:)
      double precision                          :: ph(4)
      double precision                          :: rvec(2),unodes(4,4),phibar(4,4)
      integer                                   :: k,l

      ! Node numbering:
      !
      !  4           3
      !   o---------o
      !   |         |
      !   |    +    |
      !   |         |
      !   o---------o
      !  1           2

      ! Extrapolating state to vertices of element
      ! Point 1
      rvec(1) = s%grid%x(i,j) - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%y(i,j) - s%grid%elem(i,j)%yc
      unodes(:,1) = s%grid%elem(i,j)%u + gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Point 2
      rvec(1) = s%grid%x(i+1,j) - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%y(i+1,j) - s%grid%elem(i,j)%yc
      unodes(:,2) = s%grid%elem(i,j)%u + gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Point 3
      rvec(1) = s%grid%x(i+1,j+1) - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%y(i+1,j+1) - s%grid%elem(i,j)%yc
      unodes(:,3) = s%grid%elem(i,j)%u + gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Point 4
      rvec(1) = s%grid%x(i,j+1) - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%y(i,j+1) - s%grid%elem(i,j)%yc
      unodes(:,4) = s%grid%elem(i,j)%u + gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Finding value of limiter at each vertex
      ! k - vertex, l - variables
      do k=1,4
        do l=1,4
          if (unodes(l,k)-s%grid%elem(i,j)%u(l).gt.0.0d0) then
            phibar(l,k) = min(1.0d0,(umax(i,j,l)-s%grid%elem(i,j)%u(l))/(unodes(l,k)-s%grid%elem(i,j)%u(l)))
          else if (unodes(l,k)-s%grid%elem(i,j)%u(l).lt.0.0d0) then
            phibar(l,k) = min(1.0d0,(umin(i,j,l)-s%grid%elem(i,j)%u(l))/(unodes(l,k)-s%grid%elem(i,j)%u(l)))
          else
            phibar(l,k) = 1.0d0
          end if
        end do
      end do

      ! Finding the final value of the limiter for each conserved variable
      do k=1,4
        ph(k) = minval(phibar(k,:))
      end do

    end function barth

end module solvers
