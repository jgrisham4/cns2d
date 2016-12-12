!===========================================================
! This module contains a simple solver for the 2D Euler
! equations using a second-order accurate finite volume
! method.
!===========================================================
module euler_solver
  use cgns
  use mesh_class, only : mesh
  use utils,      only : u_to_w, w_to_u
  use limiters,   only : minmod, vanleer, barth
  use flux,       only : fluxax, fluxay
  use grad,       only : compute_gradient
  use riemann,    only : roe
  implicit none
  private
  public :: solver, initialize, solve, update, write_results, write_tec

  !---------------------------------------------------------
  ! Class for solver
  !---------------------------------------------------------
  type solver
    integer                       :: ntsteps            ! Number of time steps
    integer                       :: nelem              ! Number of elements
    double precision              :: dt                 ! Time step
    double precision              :: tfinal             ! Final time
    double precision              :: g                  ! Ratio of specific heats
    double precision              :: winfty(4)          ! Freestream primitive variables
    type(mesh)                    :: grid               ! Mesh object
    double precision, allocatable :: w_history(:,:,:,:) ! primitive solution hist
  end type solver

  contains

    !---------------------------------------------------------
    ! Subroutine which initializes solution, i.e., allocates
    ! memory and sets up some important variables
    !---------------------------------------------------------
    subroutine initialize(this,m,delta_t,t_final,gam,w0,winf)
      implicit none
      type(solver),       intent(inout)           :: this
      type(mesh),         intent(in)              :: m
      double precision,   intent(in)              :: delta_t,gam,t_final
      double precision,   intent(in), allocatable :: w0(:,:,:)
      double precision,   intent(in)              :: winf(4)
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
      write (*,'(a,i5)') "number of time steps: ", this%ntsteps

      ! Allocating memory for solution history
      allocate(this%w_history(this%grid%nelemi,this%grid%nelemj,this%ntsteps+1,4),stat=allocate_err)
      if (allocate_err.ne.0) then
        print *, "Error: Can't allocate memory for this%w_history."
        stop
      end if

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
    ! Subroutine for solving to some final time
    !---------------------------------------------------------
    subroutine solve(this)
      implicit none
      type(solver), intent(inout) :: this
      integer :: i,j,k

      ! Marching in time
      do k=1,this%ntsteps

        ! Printing some information
        write(*,'(a,i4,a,es12.5)') "timestep: ", k, " t = ", dble(k-1)*this%dt

        ! Storing current solution
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            this%w_history(i,j,k,:) = this%grid%elem(i,j)%w
          end do
        end do

        ! Updating solution
        print *, "Updating..."
        call update(this)

      end do

      ! Storing last solution
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%w_history(i,j,this%ntsteps+1,:) = this%grid%elem(i,j)%w
        end do
      end do

    end subroutine solve

    !---------------------------------------------------------
    ! Subroutine for advancing the solution one step in time
    !---------------------------------------------------------
    subroutine update(this)
      implicit none
      type(solver), intent(inout)    :: this
      double precision               :: duL(4),duR(4)
      double precision, dimension(4) :: fx,fy,uextrap,wextrap
      double precision               :: umax,umin
      integer                        :: i,j,k,err
      double precision, dimension(2) :: rL,rR,r
      double precision, allocatable  :: u(:,:,:)
      double precision, allocatable  :: gradU(:,:,:)

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
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j-1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i+1,j)%u(k),this%grid%elem(i,j-1)%u(k)))
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
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k)))
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
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j+1)%u(k)))
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
                umax = max(this%grid%elem(i,j)%u(k),max(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k)))
                umin = min(this%grid%elem(i,j)%u(k),min(this%grid%elem(i-1,j)%u(k),this%grid%elem(i,j-1)%u(k)))
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
      ! Vertical faces
      do j=1,this%grid%nelemj
        do i=2,this%grid%nelemi
          this%grid%edges_v(i,j)%flux = roe(this%grid%edges_v(i,j)%uL,this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2))
        end do
      end do
      ! Horizontal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%edges_h(i,j)%flux = roe(this%grid%edges_h(i,j)%uL,this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3))
        end do
      end do

      ! Advance one step in time (forward Euler for now)
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%elem(i,j)%u = this%grid%elem(i,j)%u - this%dt/(this%grid%elem(i,j)%area)* &
            (this%grid%edges_h(i,j)%flux*this%grid%edges_h(i,j)%length + &
            this%grid%edges_h(i,j+1)%flux*this%grid%edges_h(i,j+1)%length + &
            this%grid%edges_v(i,j)%flux*this%grid%edges_v(i,j)%length + &
            this%grid%edges_v(i+1,j)%flux*this%grid%edges_v(i+1,j)%length)
        end do
      end do

    end subroutine update

    !---------------------------------------------------------
    ! Subroutine for writing results out to a file
    !---------------------------------------------------------
    subroutine write_results(this, file_name)
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
      zonename   = "Euler solution"

      ! Opening CGNS file
      call cg_open_f(file_name,CG_MODE_WRITE,idx_file,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

      ! Creating zone

      ! Closing CGNS file
      call cg_close_f(idx_file,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

    end subroutine write_results

    !---------------------------------------------------------
    ! Subroutine for writing results out to a Tecplot file
    !---------------------------------------------------------
    subroutine write_tec(this, file_name)
      implicit none
      type(solver), intent(in)  :: this
      character (len=*),intent(in) :: file_name
      integer ::i,j,k

      ! Opening Tecplot file
      open(2,file=file_name)

      ! Writing header
      write(2,'(a)') 'title="Unsteady results"'
      write(2,'(a)') 'variables="x","y","rho","u","v","E"'

      ! Writing data for each time step
      do k=1,this%ntsteps

        ! Writing header information
        write(2,'(a,i5,a,i5)') 'zone i=', this%grid%imax, ' j=', this%grid%jmax
        write(2,'(a)') 'datapacking=block'
        write(2,'(a)') 'varlocation=([3,4,5,6]=cellcentered)'
        write(2,'(a,es25.10)') 'solutiontime=', dble(k-1)*this%dt

        ! Writing x-coordinates of nodes
        do j=1,this%grid%jmax
          do i=1,this%grid%imax
            write(2,'(es25.10)',advance='no') this%grid%x(i,j)
            if (mod((j-1)*this%grid%imax+i,10).eq.0) then
              write(2,'(a)') " "
            end if
          end do
        end do
        if (mod((j-1)*this%grid%imax+i,10).ne.0) then
          write(2,'(a)') " "
        end if

        ! Writing y-coordinates of nodes
        do j=1,this%grid%jmax
          do i=1,this%grid%imax
            write(2,'(es25.10)',advance='no') this%grid%y(i,j)
            if (mod((j-1)*this%grid%imax+i,10).eq.0) then
              write(2,'(a)') " "
            end if
          end do
        end do
        if (mod((j-1)*this%grid%imax+i,10).ne.0) then
          write(2,'(a)') " "
        end if

        ! Writing density at cell-centers
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            write(2,'(es25.10)',advance='no') this%w_history(i,j,k,1)
            if (mod((j-1)*this%grid%nelemi+i,10).eq.0) then
              write(2,'(a)') " "
            end if
          end do
        end do
        if (mod((j-1)*this%grid%nelemi+i,10).ne.0) then
          write(2,'(a)') " "
        end if

        ! Writing x-velocity at cell-centers
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            write(2,'(es25.10)',advance='no') this%w_history(i,j,k,2)
            if (mod((j-1)*this%grid%nelemi+i,10).eq.0) then
              write(2,'(a)') " "
            end if
          end do
        end do
        if (mod((j-1)*this%grid%nelemi+i,10).ne.0) then
          write(2,'(a)') " "
        end if

        ! Writing y-velocity at cell-centers
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            write(2,'(es25.10)',advance='no') this%w_history(i,j,k,3)
            if (mod((j-1)*this%grid%nelemi+i,10).eq.0) then
              write(2,'(a)') " "
            end if
          end do
        end do
        if (mod((j-1)*this%grid%nelemi+i,10).ne.0) then
          write(2,'(a)') " "
        end if

        ! Writing pressure at cell-centers
        do j=1,this%grid%nelemj
          do i=1,this%grid%nelemi
            write(2,'(es25.10)',advance='no') this%w_history(i,j,k,4)
            if (mod((j-1)*this%grid%nelemi+i,10).eq.0) then
              write(2,'(a)') " "
            end if
          end do
        end do
        if (mod((j-1)*this%grid%nelemi+i,10).ne.0) then
          write(2,'(a)') " "
        end if

      end do

      ! Closing file
      close(2)

    end subroutine write_tec

end module euler_solver
