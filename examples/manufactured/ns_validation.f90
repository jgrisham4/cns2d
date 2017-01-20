module validation
  use mesh_class,   only : mesh,read_from_file,preprocess
  use solver_class, only : solver,initialize
  use temporal,     only : solve_steady
  use mms,          only : rho_e,u_e,v_e,et_e
  use utils,        only : u_to_w
  implicit none
  private
  public :: compute_error

  contains

    subroutine compute_error(dim_i,dim_j,cfl_num, dx, error_norm)
      implicit none
      integer, intent(in)             :: dim_i,dim_j
      double precision, intent(in)    :: cfl_num
      double precision, intent(inout) :: dx,error_norm
      double precision, parameter     :: g = 1.4d0
      double precision, parameter     :: R = 1.0d0
      double precision                :: dt,t_final,xc,yc,xerr,yerr
      double precision                :: rhoc,uc,vc,etc,utmp(4),winfty(4),tol,cfl
      double precision, allocatable   :: w0(:,:,:)
      integer                         :: aerr,i,j,bcs(4),niter,niterfo
      character(len=30)               :: ext,gridname

      ! Creating mesh and solver objects
      type(mesh)   :: grid
      type(solver) :: ns_solver

      write(gridname,'(a,i0,a,i0,a)') "grid_", dim_i, "x", dim_j, ".cgns"

      ! Reading mesh
      call read_from_file(grid,gridname)

      ! Preprocessing grid
      call preprocess(grid)
      dx = sqrt(grid%elem(1,1)%area)

      ! Allocating memory for the initial guess
      allocate(w0(grid%nelemi,grid%nelemj,4),stat=aerr)
      if (aerr.ne.0) then
        print *, "Can't allocate memory for w0."
        stop
      end if

      ! Setting time step and final time
      dt      = 1.0d0
      t_final = 0.0d0
      cfl     = cfl_num
      tol     = 1.0e-10
      niterfo = 200   ! First-order iterations
      niter   = 1000

      ! Setting initial guess using random noise added to the initial guess
      winfty(:) = 0.0d0
      do j=1,grid%nelemj
        do i=1,grid%nelemi

          ! Generating pseudo-random noise which will be added to coordinates of cell centers
          call random_number(xerr)  ! Returns number 0 <= xerr < 1
          call random_number(yerr)
          xerr = xerr*0.02d0 - 0.01d0
          yerr = yerr*0.02d0 - 0.01d0
          xerr = 0.0d0
          yerr = 0.0d0

          ! Computing the exact solution at the cell center (with noise added)
          xc   = grid%elem(i,j)%xc + xerr
          yc   = grid%elem(i,j)%yc + yerr
          rhoc = rho_e(xc,yc)
          uc   = u_e(xc,yc)
          vc   = v_e(xc,yc)
          etc  = et_e(xc,yc)

          ! Constructing the vector of conserved variables
          utmp(1) = rhoc
          utmp(2) = rhoc*uc
          utmp(3) = rhoc*vc
          utmp(4) = rhoc*etc

          ! Setting primitive state
          w0(i,j,:) = u_to_w(utmp,g)

        end do
      end do

      ! Setting boundary conditions
      bcs(1) = 2000
      bcs(2) = 2000
      bcs(3) = 2000
      bcs(4) = 2000

      ! Initializing solver
      call initialize(ns_solver,grid,dt,t_final,g,R,w0,winfty,bcs,"none",.true.,niter,niterfo,tol,cfl)

      ! Solving the problem
      call solve_steady(ns_solver)

      ! Computing the L2 norm of the error
      error_norm = 0.0d0
      do j=1,ns_solver%grid%nelemj
        do i=1,ns_solver%grid%nelemi
          error_norm = error_norm + (ns_solver%grid%elem(i,j)%u(1) - &
            rho_e(ns_solver%grid%elem(i,j)%xc,ns_solver%grid%elem(i,j)%yc))**2*&
            ns_solver%grid%elem(i,j)%area
        end do
      end do
      error_norm = sqrt(error_norm)

    end subroutine compute_error

end module validation

program ns_validation
  use validation
  implicit none
  double precision, dimension(3) :: err,nelem,xv,yv,dx
  double precision               :: A(3,2),xstar(2),b(3)
  double precision               :: Am(2,2),bm(2)

  ! Computing the error in each solution (only looking at density for now)
  !call compute_error( 41,21,0.001d0,err(1),dx(1))
  !call compute_error( 81,41, 0.01d0,err(2),dx(2))
  call compute_error(161,81, 0.5d0,err(3),dx(3))
  !nelem(1) = 1.0d0/sqrt(40.0d0*20.0d0)
  !nelem(2) = 1.0d0/sqrt(80.0d0*40.0d0)
  !nelem(3) = 1.0d0/sqrt(160.0d0*80.0d0)

  ! Finding the slope
  !xv = log10(dx)
  !yv = log10(err)
  !A(1,1) = xv(1)
  !A(1,2) = 1.0d0
  !A(2,1) = xv(2)
  !A(2,2) = 1.0d0
  !A(3,1) = xv(3)
  !A(3,2) = 1.0d0
  !b = yv
  !Am = matmul(transpose(A),A)
  !bm = matmul(transpose(A),b)
  !xstar(1) = (bm(1)*Am(2,2) - Am(1,2)*bm(2))/(Am(1,1)*Am(2,2) - Am(1,2)*Am(2,1))
  !print *, "order of accuracy for continuity is ", xstar(1)

end program ns_validation
