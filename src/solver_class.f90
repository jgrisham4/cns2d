!===============================================================================
! This module contains the data structure for the solver.  It contains field
! data in the mesh member and solver settings.  There are methods for
! initializing the solver and for writing results to a file.
!
! Author: James Grisham
! Date: 01-01-2017
!===============================================================================

module solver_class
  use cgns
  use mesh_class,   only : mesh
  use utils,        only : w_to_u
  implicit none
  public  :: solver,initialize,write_results_cgns,write_results_tec

  !-----------------------------------------------------------------------------
  ! Class for solver
  !-----------------------------------------------------------------------------
  type solver
    integer               :: ntsteps   ! Number of time steps
    integer               :: niter     ! Number of iterations used for steady flow solves
    integer               :: niterfo   ! Number of first-order iterations
    integer, dimension(4) :: bcids     ! BC identifiers
    double precision      :: dt        ! Time step
    double precision      :: cfl       ! cfl number - only used for steady flows
    double precision      :: tol       ! Tolerance used to monitor convergence to steady state
    double precision      :: tfinal    ! Final time
    double precision      :: g         ! Ratio of specific heats
    double precision      :: R         ! Ideal gas constant for air
    double precision      :: winfty(4) ! Freestream primitive variables
    logical               :: is_visc   ! Boolean variable to turn on viscous terms
    character (len=30)    :: limiter   ! Name of slope limiter ("none" or "barth")
    type(mesh)            :: grid      ! Mesh object
  end type solver

  contains

    !---------------------------------------------------------------------------
    ! Subroutine which initializes solution, i.e., allocates
    ! memory and sets up some important variables
    !---------------------------------------------------------------------------
    subroutine initialize(this,m,delta_t,t_final,gam,R,w0,winf,bcidents,lim,visc,niter,nfo,tol,cfl)
      implicit none
      type(solver),     intent(inout)           :: this
      type(mesh),       intent(in)              :: m
      double precision, intent(in)              :: delta_t,gam,t_final,R
      double precision, intent(in), allocatable :: w0(:,:,:)
      double precision, intent(in)              :: winf(4)
      integer,          intent(in)              :: bcidents(4)
      character (len=*),intent(in)              :: lim
      logical,          intent(in)              :: visc
      integer,          intent(in)              :: niter
      integer,          intent(in)              :: nfo
      double precision, intent(in)              :: tol  ! Tolerance used to monitor convergence
      double precision, intent(in)              :: cfl
      double precision, allocatable             :: u0(:,:,:)
      integer                                   :: allocate_err,i,j

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
      this%is_visc = visc
      this%niter   = niter
      this%niterfo = nfo
      this%tol     = tol
      this%cfl     = cfl
      !write (*,'(a,i7)') "number of time steps: ", this%ntsteps

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


    !---------------------------------------------------------------------------
    ! Subroutine for writing results out to a CGNS file
    !---------------------------------------------------------------------------
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

    !---------------------------------------------------------------------------
    ! Subroutine for writing results out to a Tecplot file
    !---------------------------------------------------------------------------
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

end module solver_class
