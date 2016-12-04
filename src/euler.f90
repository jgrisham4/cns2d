!===========================================================
! This module contains a simple solver for the 2D Euler 
! equations using a second-order accurate finite volume 
! method.
!===========================================================
module euler
  use class_mesh, only : mesh
  use utils,      only : u_to_w, w_to_u
  use limiters,   only : minmod, vanleer, venkatakrishnan
  use flux,       only : inviscid_flux
  use riemann,    only : Roe
  implicit none
  private
  public :: elem_data, solver, initialize, solve, update, write_results

  !*********************************************************
  ! Class for solver
  !*********************************************************
  type solver
    integer                       :: ntsteps            ! Number of time steps
    integer                       :: nelem              ! Number of elements
    double precision              :: dt                 ! Time step
    double precision              :: g                  ! Ratio of specific heats
    type(mesh)                    :: grid               ! Mesh object
    type(elem_data), allocatable  :: elem(:)            ! Element data
    double precision, allocatable :: w_history(:,:,:,:) ! w solution hist
  end type solver

  contains 

    !---------------------------------------------------------
    ! Subroutine which initializes solution, i.e., allocates 
    ! memory and sets up some important variables
    ! 
    ! NOTE: The initial guess should include the values at
    ! the ghost cells... I think
    !---------------------------------------------------------
    subroutine initialize(this,m,nt,delta_t,gam,w0) 
      implicit none
      type(solver),     intent(inout)           :: this
      type(mesh),       intent(in)              :: m
      integer,          intent(in)              :: nt
      double precision, intent(in)              :: delta_t,gam
      double precision, intent(in), allocatable :: w0(:,:,:)
      double precision, allocatable             :: u0(:,:,:)
      integer                                   :: allocate_err,i

      ! Assigning members of the solver class
      this%grid    = m
      this%ntsteps = nt
      this%dt      = delta_t
      this%g       = gam

      ! Allocating memory for solution history
      allocate(this%w_history(this%grid%nelemi,this%grid%nelemj,nt+1,4),stat=allocate_err)
      if (allocate_err.ne.0) then
        print *, "Error: Can't allocate memory for this%w_history."
        stop
      end if

      ! Converting initial condition to conserved variables 
      ! (includes ghost cells)
      allocate(u0,mold=w0)
      do j=-1,this%grid%nelemj+2
        do i=-1,this%grid%nelemi+2
          u0(i,j,:) = w_to_u(w0(i,j,:),gam)
        end do
      end do

      ! Filling data into elem structures
      do j=-1,this%grid%nelemj+2
        do i=-1,this%nelem+2
          this%elem(i,j)%u  = u0(i,j,:)
          this%elem(i,j)%w  = w0(i,j,:)
        end do
      end do

      print *, "Done initializing solver."

    end subroutine initialize

    !---------------------------------------------------------
    ! Subroutine for solving to some final time
    !---------------------------------------------------------
    subroutine solve(this)
      implicit none
      type(solver), intent(inout) :: this
      integer :: i,j

      ! Marching in time
      do i=1,this%ntsteps

        ! Printing some information
        !write(*,'(es10.5)') dble(i-1)*this%dt
        print *, "t=", dble(i)*this%dt

        ! Storing current solution
        do j=1,this%grid%num_elements
          this%w_history(j,i,:) = this%elem(j)%w
        end do

        ! Updating solution
        call update(this)

      end do

      ! Storing current solution
      !do j=1,this%grid%num_elements
      !  this%w_history(j,i+1,:) = this%elem(j)%w
      !end do

    end subroutine solve 
    
    !---------------------------------------------------------
    ! Subroutine for advancing the solution one step in time
    !---------------------------------------------------------
    subroutine update(this)
      implicit none
      type(solver), intent(inout)   :: this
      double precision              :: dwL(3),dwR(3)
      double precision              :: wL(3),wR(3)
      double precision              :: delta(3)
      double precision              :: r(3)
      integer                       :: i,j,k,err
      double precision, allocatable :: u(:,:,:)
      double precision, allocatable :: gradU(:,:,:)
      !double precision,allocatable :: f_interface(:,:)

      ! Allocating memory for gradient of state at cell centers
      allocate(gradU(this%grid%nelemi,this%grid%nelemj,8),stat=err)
      if (err.ne.0) then
        print *, "Error: Can't allocate memory for gradU in update."
        stop
      end if
      allocate(u(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: Can't allocate memory for u in update."
        stop
      end if 

      ! Computing gradients of the state in each cell for use 
      ! in piecewise linear reconstruction

      ! Reconstructing states on left and right sides of each 
      ! interior cell (CORRECT)
      ! UNIFORM ONLY
      do i=1,this%nelem-1
        
        ! Finding slopes on left and right
        ! UNIFORM
        dwL = this%elem(i)%w - this%elem(i-1)%w
        dwR = this%elem(i+1)%w - this%elem(i)%w
        r = dwR/dwL

        ! Reconstructing primitive states on left and right of element
        do j=1,3
          this%elem(i)%wL(j) = this%elem(i)%w(j) + & 
            0.5d0*(this%elem(i)%w(j)-this%elem(i-1)%w(j))*minmod(r(j))
          this%elem(i)%wR(j) = this%elem(i)%w(j) - &
            0.5d0*(this%elem(i+1)%w(j)-this%elem(i)%w(j))*minmod(1.0d0/r(j))
        end do


      end do

      ! The current boundary condition is no gradient at either end.
      ! This is enforced by copying the state from the interior of the domain
      ! to the ghost cells.  
      this%elem(0)%wL = this%elem(1)%wR
      this%elem(this%nelem-1)%wR = this%elem(this%nelem-2)%wL
      this%elem(0)%w = this%elem(1)%w
      !this%elem(0)%w(2) = -this%elem(1)%w(2)   ! v_ghost = -v_1
      this%elem(this%nelem-1)%w = this%elem(this%nelem-2)%w
      !this%elem(this%nelem)%w(2) = -this%elem(this%nelem-1)%w(2)

      ! Must iterate through all the interior interfaces and solve 
      ! the Riemann problem to find the fluxes
      do i=1,this%nelem-1
        f_interface(:,i) = roe_flux2(this%elem(i-1)%wL,this%elem(i)%wR,this%g)
        !f_interface(:,i) = roe_flux(this%elem(i-1)%wL,this%elem(i)%wR,this%g)
        !print *, f_interface(2,i)
      end do

      ! Advance one step in time (forward Euler for now)
      do i=1,this%nelem-2
        this%elem(i)%u0 = this%elem(i)%u
        this%elem(i)%u = this%elem(i)%u - this%dt/(this%grid%dx)* &
          (f_interface(:,i+1) - f_interface(:,i))
        this%elem(i)%w = u_to_w(this%elem(i)%u,this%g)
      end do

    end subroutine update

    !---------------------------------------------------------
    ! Subroutine for writing results out to a file
    !---------------------------------------------------------
    subroutine write_results(this, file_name)
      implicit none
      type(solver), intent(inout) :: this
      character (len=*)           :: file_name
      integer                     :: i,j

      ! Opening file
      open(3,file="inputs.dat")
      write(3,*) this%grid%num_elements, this%ntsteps, this%dt
      
      ! Writing cell centers
      do i=1,this%grid%num_elements
        write(3,*) this%grid%xc(i)
      end do
      
      ! Closing file
      close(3)

      ! Writing density
      open(3,file="rho.dat")
      do j=1,this%ntsteps
        do i=1,this%grid%num_elements
          write(3,*) this%w_history(i,j,1)
        end do
      end do
      close(3)
        
      ! Writing velocity
      open(3,file="u.dat")
      do j=1,this%ntsteps
        do i=1,this%grid%num_elements
          write(3,*) this%w_history(i,j,2)
        end do
      end do
      close(3)
        
      ! Writing pressure
      open(3,file="p.dat")
      do j=1,this%ntsteps
        do i=1,this%grid%num_elements
          write(3,*) this%w_history(i,j,3)
        end do
      end do
      close(3)


    end subroutine write_results

end module class_solver
