!===============================================================================
! This module contains the subroutines which are used to march the solver
! in time.  All explicit for now.
!
! Author: James Grisham
! Date: 01/13/2017
!===============================================================================

module temporal
  use mesh_class,   only : compute_max_timesteps_inv,compute_max_timesteps_visc
  use solver_class, only : solver,write_results_tec
  use euler,        only : residual_inv_fo, residual_inv
  use navierstokes, only : residual_visc_fo, residual_visc
  use linalg,       only : norml2
  use acceleration, only : irs_upwind

  ! Subroutines and functions
  !public :: solve_feuler, solve_rk4, solve_steady, solve_steady_fe
  public :: solve_feuler,solve_rk4,solve_steady

  contains

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
      if (this%is_visc.eq..true.) then
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
          call residual_visc(this,resid)

          ! Advancing in time
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt*resid(i,j,:)
            end do
          end do

        end do
      else
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
      end if

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
      if (this%is_visc.eq..true.) then
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
          call residual_visc(this,resid)
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt/4.0d0*resid(i,j,:)
            end do
          end do

          ! Stage 2
          call residual_visc(this,resid)
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt/3.0d0*resid(i,j,:)
            end do
          end do

          ! Stage 3
          call residual_visc(this,resid)
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt/2.0d0*resid(i,j,:)
            end do
          end do

          ! Stage 4
          call residual_visc(this,resid)
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt*resid(i,j,:)
            end do
          end do

        end do
      else
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
      end if

    end subroutine solve_rk4

    !---------------------------------------------------------
    ! Subroutine for solving a steady problem -- a 3 stage
    ! Runge-Kutta scheme is used to advance in time.
    !---------------------------------------------------------
    subroutine solve_steady(this)
      implicit none
      type(solver), intent(inout)    :: this
      double precision, allocatable  :: resid(:,:,:),eqn_resids(:,:)
      double precision, dimension(4) :: eqn_resids0
      double precision, allocatable :: rbar(:,:,:)
      character (len=30)             :: tecname
      integer                        :: i,j,k,l,aer

      ! Parameter used in upwind implicit residual smoothing
      double precision, parameter :: epsln = 0.10d0

      ! Parameters used in multi-stage time-stepping
      !double precision, parameter   :: a1 = 0.1918d0
      !double precision, parameter   :: a2 = 0.4929d0   ! These three are for 2nd-order accuracy
      !double precision, parameter   :: a3 = 1.0000d0
      double precision, parameter   :: a1 = 0.1481d0
      double precision, parameter   :: a2 = 0.4000d0
      double precision, parameter   :: a3 = 1.0000d0

      ! Allocating memory for the solution and the residual
      allocate(resid(this%grid%nelemi,this%grid%nelemj,4),stat=aer)
      if (aer.ne.0) then
        print *, "Error: can't allocate memory for resid in solve_steady."
        stop
      end if
      allocate(eqn_resids(this%niter+this%niterfo,4),stat=aer)
      if (aer.ne.0) then
        print *, "Error: can't allocate memory for eqn_resids in solve_steady."
        stop
      end if
      allocate(rbar,mold=resid,stat=aer)
      if (aer.ne.0) then
        print *, "Error: can't allocate memory for rbar in solve_steady."
        stop
      end if

      ! Writing initial solution
      tecname = "initial.tec"
      print *, "Writing data to ", tecname
      call write_results_tec(this,tecname)

      ! Computing time steps
      if (this%is_visc.eq..true.) then
        call compute_max_timesteps_visc(this%grid,this%g,this%R,this%cfl)
      else
        call compute_max_timesteps_inv(this%grid,this%g,this%R,this%cfl)
      end if

      ! Setting initial residuals
      eqn_resids(1,:) = 1.0d0
      eqn_resids0(:) = 1.0d0
      !write(*,'(a,i6)',advance='no') "iteration: ", 0
      !write(*,'(a,4(es12.5,x))') " residuals: ", eqn_resids(1,:)

      ! Marching in time
      k = 1
      if (this%is_visc.eq..true.) then

        ! First-order iterations
        do while (k.le.this%niterfo)

          ! Copying old solution
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              this%grid%elem(i,j)%u0 = this%grid%elem(i,j)%u
            end do
          end do

          ! Stage 1
          call residual_visc_fo(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a1*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a1*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Updating local time steps
          call compute_max_timesteps_visc(this%grid,this%g,this%R,this%cfl)

          ! Stage 2
          call residual_visc_fo(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a2*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a2*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Updating local time steps
          call compute_max_timesteps_visc(this%grid,this%g,this%R,this%cfl)

          ! Stage 3
          call residual_visc_fo(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a3*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a3*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Computing the residuals
          do l=1,4
            eqn_resids(k,l) = norml2(this%grid,resid(:,:,l))
          end do

          ! Getting initial residual
          !if (k.eq.1) then
          !  eqn_resids0(:) = eqn_resids(k,:)
          !end if

          ! Updating local time steps
          call compute_max_timesteps_visc(this%grid,this%g,this%R,this%cfl)

          ! Printing some information
          write(*,'(a,i6)',advance='no') "iteration: ", k
          write(*,'(a,4(es12.5,x))') " residuals: ", eqn_resids(k,:)/eqn_resids0

          ! Incrementing counter
          k = k + 1

        end do

        ! Switching to second-order spatial accuracy
        do while (k.le.(this%niter+this%niterfo))
        !do while ((any(eqn_resids(k,:).ge.this%tol)).and.(k.le.this%niter))

          ! Copying old solution
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              this%grid%elem(i,j)%u0 = this%grid%elem(i,j)%u
            end do
          end do

          ! Stage 1
          call residual_visc(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a1*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a1*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Updating local time steps
          call compute_max_timesteps_visc(this%grid,this%g,this%R,this%cfl)

          ! Stage 2
          call residual_visc(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a2*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a2*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Updating local time steps
          call compute_max_timesteps_visc(this%grid,this%g,this%R,this%cfl)

          ! Stage 3
          call residual_visc(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a3*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a3*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Computing the residuals
          do l=1,4
            eqn_resids(k,l) = norml2(this%grid,resid(:,:,l))
          end do

          ! Getting initial residual
          !if (k.eq.1) then
          !  eqn_resids0(:) = eqn_resids(k,:)
          !end if

          ! Updating local time steps
          call compute_max_timesteps_visc(this%grid,this%g,this%R,this%cfl)

          ! Printing some information
          write(*,'(a,i6)',advance='no') "iteration: ", k
          write(*,'(a,4(es12.5,x))') " residuals: ", eqn_resids(k,:)/eqn_resids0

          ! Incrementing counter
          k = k + 1

        end do
      else

        ! Inviscid first-order spatial accuracy
        do while (k.le.this%niterfo)

          ! Copying old solution
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              this%grid%elem(i,j)%u0 = this%grid%elem(i,j)%u
            end do
          end do

          ! Stage 1
          call residual_inv_fo(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a1*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a1*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Updating local time steps
          call compute_max_timesteps_inv(this%grid,this%g,this%R,this%cfl)

          ! Stage 2
          call residual_inv_fo(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a2*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a2*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Updating local time steps
          call compute_max_timesteps_inv(this%grid,this%g,this%R,this%cfl)

          ! Stage 3
          call residual_inv_fo(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a3*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a3*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Computing the residuals
          do l=1,4
            eqn_resids(k,l) = norml2(this%grid,resid(:,:,l))
          end do

          ! Getting initial residual
          !if (k.eq.1) then
          !  eqn_resids0(:) = eqn_resids(k,:)
          !end if

          ! Updating local time steps
          call compute_max_timesteps_inv(this%grid,this%g,this%R,this%cfl)

          ! Printing some information
          write(*,'(a,i6)',advance='no') "iteration: ", k
          write(*,'(a,4(es12.5,x))') " residuals: ", eqn_resids(k,:)/eqn_resids0

          ! Incrementing counter
          k = k + 1

        end do

        ! Switching to second-order spatial accuracy
        !do while ((all(eqn_resids(k,:).ge.this%tol)).and.(k.le.this%niter))
        do while (k.le.(this%niter+this%niterfo))

          ! Copying old solution
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              this%grid%elem(i,j)%u0 = this%grid%elem(i,j)%u
            end do
          end do

          ! Stage 1
          call residual_inv(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a1*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a1*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Updating local time steps
          call compute_max_timesteps_inv(this%grid,this%g,this%R,this%cfl)

          ! Stage 2
          call residual_inv(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a2*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a2*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Updating local time steps
          call compute_max_timesteps_inv(this%grid,this%g,this%R,this%cfl)

          ! Stage 3
          call residual_inv(this,resid)
          call irs_upwind(this%grid,epsln,resid,this%g,rbar)
          resid = rbar
          do j=1,this%grid%nelemj
            do i=1,this%grid%nelemi
              !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a3*this%dt*resid(i,j,:)
              this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + a3*this%grid%elem(i,j)%dt_max*resid(i,j,:)
            end do
          end do

          ! Computing the residuals
          do l=1,4
            eqn_resids(k,l) = norml2(this%grid,resid(:,:,l))
          end do

          ! Getting initial residual
          !if (k.eq.1) then
          !  eqn_resids0(:) = eqn_resids(k,:)
          !end if

          ! Updating local time steps
          call compute_max_timesteps_inv(this%grid,this%g,this%R,this%cfl)

          ! Printing some information
          write(*,'(a,i6)',advance='no') "iteration: ", k
          write(*,'(a,4(es12.5,x))') " residuals: ", eqn_resids(k,:)/eqn_resids0

          ! Incrementing counter
          k = k + 1

        end do
      end if

      ! Writing final results to file
      tecname = "solution.tec"
      print *, "Writing final solution to ", tecname
      call write_results_tec(this,tecname)

      ! Writing convergence history to file
      open(2,file="residuals.dat")
      do i=1,(k-1)
        write(2,'(i6,4(es12.5,x))') i, eqn_resids(i,:)/eqn_resids0
      end do
      close(2)

    end subroutine solve_steady

    !---------------------------------------------------------
    ! Subroutine for solving a steady problem -- the forward
    ! Euler method is used to advance in time.
    !---------------------------------------------------------
    !subroutine solve_steady_fe(this)
    !  implicit none
    !  type(solver), intent(inout)    :: this
    !  double precision, allocatable  :: resid(:,:,:),eqn_resids(:,:)
    !  double precision, dimension(4) :: eqn_resids0
    !  !double precision, allocatable :: rbar(:,:,:)
    !  character (len=30)             :: tecname
    !  integer                        :: i,j,k,l,aer

    !  ! Parameter used in upwind implicit residual smoothing
    !  double precision, parameter :: epsln = 4.0d0

    !  ! Parameters used in multi-stage time-stepping
    !  double precision, parameter   :: a1 = 0.1481d0
    !  double precision, parameter   :: a2 = 0.4000d0
    !  double precision, parameter   :: a3 = 1.0000d0

    !  ! Allocating memory for the solution and the residual
    !  allocate(resid(this%grid%nelemi,this%grid%nelemj,4),stat=aer)
    !  if (aer.ne.0) then
    !    print *, "Error: can't allocate memory for resid in solve_steady."
    !    stop
    !  end if
    !  allocate(eqn_resids(this%niter+this%niterfo,4),stat=aer)
    !  if (aer.ne.0) then
    !    print *, "Error: can't allocate memory for eqn_resids in solve_steady."
    !    stop
    !  end if
    !  !allocate(rbar,mold=resid,stat=aer)
    !  !if (aer.ne.0) then
    !  !  print *, "Error: can't allocate memory for rbar in solve_steady."
    !  !  stop
    !  !end if

    !  ! Writing initial solution
    !  tecname = "initial.tec"
    !  print *, "Writing data to ", tecname
    !  call write_results_tec(this,tecname)

    !  ! Computing time steps
    !  !if (this%is_visc.eq..true.) then
    !  !  call compute_max_timesteps_visc(this%grid,this%g,this%R,this%cfl)
    !  !else
    !  !  call compute_max_timesteps_inv(this%grid,this%g,this%R,this%cfl)
    !  !end if

    !  ! Setting initial residuals
    !  eqn_resids(1,:) = 1.0d0
    !  eqn_resids0(:) = 1.0d0
    !  write(*,'(a,i6)',advance='no') "iteration: ", 0
    !  write(*,'(a,4(es12.5,x))') " residuals: ", eqn_resids(1,:)/eqn_resids0

    !  ! Marching in time
    !  k = 1
    !  if (this%is_visc.eq..true.) then

    !    ! First-order iterations
    !    do while (k.le.this%niterfo)

    !      ! Copying old solution
    !      do j=1,this%grid%nelemj
    !        do i=1,this%grid%nelemi
    !          this%grid%elem(i,j)%u0 = this%grid%elem(i,j)%u
    !        end do
    !      end do

    !      ! Advancing the solution in time
    !      call residual_visc_fo(this,resid)
!   !       call irs_upwind(this%grid,epsln,resid,this%g,rbar)
!   !       resid = rbar
    !      do j=1,this%grid%nelemj
    !        do i=1,this%grid%nelemi
    !          !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%grid%elem(i,j)%dt_max*resid(i,j,:)
    !          this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt*resid(i,j,:)
    !        end do
    !      end do

    !      ! Computing the residuals
    !      do l=1,4
    !        eqn_resids(k,l) = norml2(this%grid,resid(:,:,l))
    !      end do

    !      ! Getting initial residual used for normalization
    !      if (k.eq.1) then
    !        eqn_resids0(:) = eqn_resids(k,:)
    !      end if

    !      ! Updating local time steps
    !      !call compute_max_timesteps_visc(this%grid,this%g,this%R,this%cfl)

    !      ! Printing some information
    !      write(*,'(a,i6)',advance='no') "iteration: ", k
    !      write(*,'(a,4(es12.5,x))') " residuals: ", eqn_resids(k,:)/eqn_resids0

    !      ! Incrementing counter
    !      k = k + 1

    !    end do

    !    ! Switching to second-order spatial accuracy
    !    do while (k.le.(this%niter+this%niterfo))
    !    !do while ((any(eqn_resids(k,:).ge.this%tol)).and.(k.le.this%niter))

    !      ! Copying old solution
    !      do j=1,this%grid%nelemj
    !        do i=1,this%grid%nelemi
    !          this%grid%elem(i,j)%u0 = this%grid%elem(i,j)%u
    !        end do
    !      end do

    !      ! Advancing in time using forward Euler
    !      call residual_visc(this,resid)
!   !       call irs_upwind(this%grid,epsln,resid,this%g,rbar)
!   !       resid = rbar
    !      do j=1,this%grid%nelemj
    !        do i=1,this%grid%nelemi
    !          !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%grid%elem(i,j)%dt_max*resid(i,j,:)
    !          this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt*resid(i,j,:)
    !        end do
    !      end do

    !      ! Computing the residuals
    !      do l=1,4
    !        eqn_resids(k,l) = norml2(this%grid,resid(:,:,l))
    !      end do

    !      ! Getting initial residual used for normalization
    !      if (k.eq.1) then
    !        eqn_resids0(:) = eqn_resids(k,:)
    !      end if

    !      ! Updating local time steps
    !      !call compute_max_timesteps_visc(this%grid,this%g,this%R,this%cfl)

    !      ! Printing some information
    !      write(*,'(a,i6)',advance='no') "iteration: ", k
    !      write(*,'(a,4(es12.5,x))') " residuals: ", eqn_resids(k,:)/eqn_resids0

    !      ! Incrementing counter
    !      k = k + 1

    !    end do
    !  else

    !    ! Inviscid first-order spatial accuracy
    !    do while (k.le.this%niterfo)

    !      ! Copying old solution
    !      do j=1,this%grid%nelemj
    !        do i=1,this%grid%nelemi
    !          this%grid%elem(i,j)%u0 = this%grid%elem(i,j)%u
    !        end do
    !      end do

    !      ! Advancing the solution in time using forward Euler
    !      call residual_inv_fo(this,resid)
!   !       call irs_upwind(this%grid,epsln,resid,this%g,rbar)
!   !       resid = rbar
    !      do j=1,this%grid%nelemj
    !        do i=1,this%grid%nelemi
    !          !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%grid%elem(i,j)%dt_max*resid(i,j,:)
    !          this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt*resid(i,j,:)
    !        end do
    !      end do

    !      ! Computing the residuals
    !      do l=1,4
    !        eqn_resids(k,l) = norml2(this%grid,resid(:,:,l))
    !      end do

    !      ! Getting initial residual used for normalization
    !      if (k.eq.1) then
    !        eqn_resids0(:) = eqn_resids(k,:)
    !      end if

    !      ! Updating local time steps
    !      !call compute_max_timesteps_inv(this%grid,this%g,this%R,this%cfl)

    !      ! Printing some information
    !      write(*,'(a,i6)',advance='no') "iteration: ", k
    !      write(*,'(a,4(es12.5,x))') " residuals: ", eqn_resids(k,:)/eqn_resids0

    !      ! Incrementing counter
    !      k = k + 1

    !    end do

    !    ! Switching to second-order spatial accuracy
    !    !do while ((all(eqn_resids(k,:).ge.this%tol)).and.(k.le.this%niter))
    !    do while (k.le.(this%niter+this%niterfo))

    !      ! Copying old solution
    !      do j=1,this%grid%nelemj
    !        do i=1,this%grid%nelemi
    !          this%grid%elem(i,j)%u0 = this%grid%elem(i,j)%u
    !        end do
    !      end do

    !      ! Advancing the solution in time using forward Euler
    !      call residual_inv(this,resid)
!   !       call irs_upwind(this%grid,epsln,resid,this%g,rbar)
!   !       resid = rbar
    !      do j=1,this%grid%nelemj
    !        do i=1,this%grid%nelemi
    !          !this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%grid%elem(i,j)%dt_max*resid(i,j,:)
    !          this%grid%elem(i,j)%u = this%grid%elem(i,j)%u0 + this%dt*resid(i,j,:)
    !        end do
    !      end do

    !      ! Computing the residuals
    !      do l=1,4
    !        eqn_resids(k,l) = norml2(this%grid,resid(:,:,l))
    !      end do

    !      ! Getting initial residual used for normalization
    !      if (k.eq.1) then
    !        eqn_resids0(:) = eqn_resids(k,:)
    !      end if

    !      ! Updating local time steps
    !      !call compute_max_timesteps_inv(this%grid,this%g,this%R,this%cfl)

    !      ! Printing some information
    !      write(*,'(a,i6)',advance='no') "iteration: ", k
    !      write(*,'(a,4(es12.5,x))') " residuals: ", eqn_resids(k,:)/eqn_resids0

    !      ! Incrementing counter
    !      k = k + 1

    !    end do
    !  end if

    !  ! Writing final results to file
    !  tecname = "solution.tec"
    !  print *, "Writing final solution to ", tecname
    !  call write_results_tec(this,tecname)

    !  ! Writing convergence history to file
    !  open(9,file="residuals.dat")
    !  do i=1,(k-1)
    !    write(9,'(i6,4(es12.5,x))') i, eqn_resids(i,:)/eqn_resids0
    !  end do
    !  close(9)
    !  print *, "Done."

    !end subroutine solve_steady_fe

end module temporal
