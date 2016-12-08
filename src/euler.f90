!===========================================================
! This module contains a simple solver for the 2D Euler
! equations using a second-order accurate finite volume
! method.
!===========================================================
module euler
  use mesh_class, only : mesh
  use utils,      only : u_to_w, w_to_u
  use limiters,   only : minmod, vanleer
  use flux,       only : inviscid_flux
  use grad,       only : compute_gradient
  use riemann,    only : roe
  !use bc_class,   only : bc,apply
  implicit none
  private
  public :: elem_data, solver, initialize, solve, update, write_results

  !---------------------------------------------------------
  ! Class for solver
  !---------------------------------------------------------
  type solver
    integer                       :: ntsteps            ! Number of time steps
    integer                       :: nelem              ! Number of elements
    double precision              :: dt                 ! Time step
    double precision              :: g                  ! Ratio of specific heats
    double precision              :: winfty(4)          ! Freestream primitive variables
    type(mesh)                    :: grid               ! Mesh object
    type(bc)                      :: bcs(4)             ! Array of bc objects
    type(elem_data), allocatable  :: elem(:)            ! Element data
    double precision, allocatable :: w_history(:,:,:,:) ! primitive solution hist
  end type solver

  contains

    !---------------------------------------------------------
    ! Subroutine which initializes solution, i.e., allocates
    ! memory and sets up some important variables
    !---------------------------------------------------------
    subroutine initialize(this,m,nt,delta_t,gam,w0,bcnames,winf)
      implicit none
      type(solver),       intent(inout)           :: this
      type(mesh),         intent(in)              :: m
      integer,            intent(in)              :: nt
      double precision,   intent(in)              :: delta_t,gam
      double precision,   intent(in), allocatable :: w0(:,:,:)
      character (len=30), intent(in)              :: bcnames(4)
      double precision, allocatable               :: u0(:,:,:)
      double precision, dimension(4)              :: winf
      integer                                     :: allocate_err,i

      ! Assigning members of the solver class
      this%grid    = m
      this%ntsteps = nt
      this%dt      = delta_t
      this%g       = gam
      this%winfty  = winf

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
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%elem(i,j)%u  = u0(i,j,:)
          this%elem(i,j)%w  = w0(i,j,:)
        end do
      end do

      ! Setting up boundary condition objects
      do i=1,4
        bcs(i)%bcname = bcnames(i)
      end do

      print *, "Done initializing solver."

    end subroutine initialize

    !---------------------------------------------------------
    ! Subroutine for solving to some final time
    !---------------------------------------------------------
    subroutine solve(this)
      implicit none
      type(solver), intent(inout) :: this
      integer :: i,j,k

      ! Marching in time
      do k=1,this%ntsteps

        ! Printing some information
        !write(*,'(es10.5)') dble(i-1)*this%dt
        print *, "t=", dble(i)*this%dt

        ! Storing current solution
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%w_history(i,j,k,:) = this%grid%elem(i,j)%w
          end do
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
      type(solver), intent(inout)    :: this
      double precision               :: duL(4),duR(4)
      double precision               :: wL(4),wR(4)
      double precision, dimension(4) :: fx,fy,uextrap,wextrap
      double precision               :: umax,umin
      integer                        :: i,j,k,err
      double precision, dimension(2) :: rL,rR,r
      double precision, allocatable  :: u(:,:,:)
      double precision, allocatable  :: gradU(:,:,:)
      !double precision,allocatable  :: f_interface(:,:)

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
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          u(i,j,:) = this%grid%elem(i,j)%u
        end do
      end do
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

            ! Reconstructing primitive states on left of interface
            do k=1,4

              ! Finding necessary inputs to the Barth-Jespersen limiter
              if (j.eq.1) then
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
              else if (j.eq.this%grid%nelemj) then
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j-1)%u(k))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j-1)%u(k))
              else
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i,j-1)%u(k),this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i,j-1)%u(k),this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
              end if

              ! Reconstruction
              this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + barth(duL(k),this%grid%elem(i,j)%u(k),umax,umin)

            end do

          else if (i.eq.(this%grid%nelemi)) then

            ! Only need to reconstruct the right state
            rR(1) = this%grid%edges_v(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_v(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding the slope on the right side of the interface
            duR = gradU(i,j,1:4)*rR(1) + gradU(i,j,5:8)*rR(2)

            ! Reconstructing primitive states on right of interface
            do k=1,4

              ! Finding necessary inputs to the Barth-Jespersen limiter
              if (j.eq.1) then
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
              else if (j.eq.this%grid%nelemj) then
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k))
              else
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i,j-1)%u(k),this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i,j-1)%u(k),this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
              end if

              ! Reconstruction
              this%grid%edges_v(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + barth(duR(k),this%grid%elem(i,j)%u(k),umax,umin)

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

              ! Finding necessary inputs to the Barth-Jespersen limiter
              umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k),this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
              umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k),this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))

              ! Reconstruction
              this%grid%edges_v(i+1,j)%uL(k) = this%grid%elem(i,j)%u(k) + barth(duL(k),this%grid%elem(i,j)%u(k),umax,umin)
              this%grid%edges_v(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + barth(duR(k),this%grid%elem(i,j)%u(k),umax,umin)

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

              ! Finding necessary inputs to the Barth-Jespersen limiter
              if (i.eq.1) then
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
              else if (i.eq.this%grid%nelemi) then
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j+1)%u(k))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j+1)%u(k))
              else
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i-1,j)%u(k),this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i-1,j)%u(k),this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
              end if

              ! Reconstruction
              this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + barth(duL(k),this%grid%elem(i,j)%u(k),umax,umin)

            end do

          else if (j.eq.(this%grid%nelemj)) then

            ! Only need to reconstruct the right state
            rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc

            ! Finding the slope on the right side of the interface
            duR = gradU(i,j,1:4)*rR(1) + gradU(i,j,5:8)*rR(2)

            ! Reconstructing primitive states on right of interface
            do k=1,4

              ! Finding necessary inputs to the Barth-Jespersen limiter
              if (i.eq.1) then
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j-1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j-1)%u(k)))
              else if (i.eq.this%grid%nelemi) then
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k))
              else
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i,i+1)%u(k),this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i,i+1)%u(k),this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k)))
              end if

              ! Reconstruction
              this%grid%edges_h(i,j)%uR(k) = this%grid%elem(i,j)%u(k) + barth(duR(k),this%grid%elem(i,j)%u(k),umax,umin)

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

              ! Finding necessary inputs to the Barth-Jespersen limiter
              umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k),this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
              umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k),this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j+1)%u(k)))

              ! Reconstruction
              this%grid%edges_h(i,j+1)%uL(k) = this%grid%elem(i,j)%u(k) + barth(duL(k),this%grid%elem(i,j)%u(k),umax,umin)
              this%grid%edges_h(i,j)%uR(k)   = this%grid%elem(i,j)%u(k) + barth(duR(k),this%grid%elem(i,j)%u(k),umax,umin)

            end do
          end if

        end do
      end do

      ! Enforcing farfield boundary conditions on the bottom and left boundaries
      ! Sides of an element are numbered as follows
      !         3
      !   o-----------o
      !   |           |
      ! 4 |           | 2
      !   |           |
      !   o-----------o
      !         1
      !
      fx = fluxax(this%winfty,this%g)
      fy = fluxay(this%winfty,this%g)
      do i=1,this%grid%nelemi
        this%grid%edges_h(i,1)%flux = fx*this%grid%elem(i,1)%n(1,1) + fy*this%grid%elem(i,1)%n(2,1)
      end do
      do j=1,this%grid%nelemj
        this%grid%edges_v(1,j)%flux = fx*this%grid%elem(1,j)%n(1,4) + fy*this%grid%elem(1,j)%n(2,4)
      end do

      ! Enforcing extrapolate bcs on the right and top boundaries
      j = this%grid%nelemj
      do i=1,this%grid%nelemi+1
        r(1) = this%grid%edges_h(i,j+1)%xm - this%grid%elem(i,j)%xc
        r(2) = this%grid%edges_h(i,j+1)%ym - this%grid%elem(i,j)%yc
        uextrap = this%grid%elem(i,j)%u + gradU(i,j,1:4)*r(1) + gradU(i,j,5:8)*r(2)
        wextrap = u_to_w(uextrap,this%g)
        fx = fluxax(wextrap,this%g)
        fy = fluxay(wextrap,this%g)
        this%grid%edges_h(i,j+1)%flux = -fx*this%grid%elem(i,j)%n(1,2) - fy*this%grid%elem(i,j)%n(2,2)
      end do
      i = this%grid%nelemi
      do j=1,this%grid%nelemj+1
        r(1) = this%grid%edges_v(i+1,j)%xm - this%grid%elem(i,j)%xc
        r(2) = this%grid%edges_v(i+1,j)%ym - this%grid%elem(i,j)%yc
        uextrap = this%grid%elem(i,j)%u + gradU(i,j,1:4)*r(1) + gradU(i,j,5:8)*r(2)
        wextrap = u_to_w(uextrap,this%g)
        fx = fluxax(wextrap,this%g)
        fy = fluxay(wextrap,this%g)
        this%grid%edges_v(i+1,j)%flux = -fx*this%grid%elem(i,j)%n(1,3) - fy*this%grid%elem(i,j)%n(2,3)
      end do

      ! Must now iterate through all the interior interfaces and solve
      ! the Riemann problem to find the fluxes
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
        end do
      end do

      ! This was in 1D
      !do i=1,this%nelem-1
      !  f_interface(:,i) = roe_flux2(this%elem(i-1)%wL,this%elem(i)%wR,this%g)
      !end do

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
