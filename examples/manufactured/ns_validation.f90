program ns_validation
  use mesh_class, only : mesh,read_from_file,preprocess
  use solvers,    only : solver,initialize,solve_steady
  use mms,        only : rho_e,u_e,v_e,et_e
  use utils,      only : u_to_w
  implicit none
  double precision, parameter   :: g = 1.4d0
  double precision, parameter   :: R = 1.0d0
  double precision              :: dt,t_final,xc,yc,xerr,yerr
  double precision              :: rhoc,uc,vc,etc,utmp(4),winfty(4),tol,cfl
  double precision, allocatable :: w0(:,:,:)
  integer                       :: aerr,i,j,bcs(4),niter

  ! Creating mesh and solver objects
  type(mesh)   :: grid
  type(solver) :: ns_solver

  ! Reading mesh
  call read_from_file(grid,"grid_41x21.cgns")
  !call read_from_file(grid,"grid_161x81.cgns")

  ! Preprocessing grid
  call preprocess(grid)

  ! Allocating memory for the initial guess
  allocate(w0(grid%nelemi,grid%nelemj,4),stat=aerr)
  if (aerr.ne.0) then
    print *, "Can't allocate memory for w0."
    stop
  end if

  ! Setting time step and final time
  cfl     = 0.0d0  ! Not used for now
  t_final = 0.0d0  ! Not used
  tol     = 1.0e-10
  dt      = 1.0e-8
  niter   = 50000

  ! Setting initial guess using random noise added to the initial guess
  winfty(:) = 0.0d0
  do j=1,grid%nelemj
    do i=1,grid%nelemi

      ! Generating pseudo-random noise which will be added to coordinates of cell centers
      call random_number(xerr)  ! Returns number 0 <= xerr < 1
      call random_number(yerr)
      xerr = xerr*0.02d0
      yerr = yerr*0.02d0

      ! Computing the exact solution at the cell center (with noise added)
      xc   = grid%elem(i,j)%xc + xerr
      yc   = grid%elem(i,j)%yc + yerr
      rhoc = rho_e(xc,yc)*1.01
      uc   = u_e(xc,yc)*1.01
      vc   = v_e(xc,yc)*1.01
      etc  = et_e(xc,yc)*1.01

      ! Constructing the vector of conserved variables
      utmp(1) = rhoc
      utmp(2) = rhoc*uc
      utmp(3) = rhoc*vc
      utmp(4) = etc
      !utmp(1) = 0.5d0
      !utmp(2) = 0.5d0*1.0d0
      !utmp(3) = 0.5d0*0.1d0
      !utmp(4) = 0.5d0*0.5d0

      ! Setting primitive state
      w0(i,j,:) = u_to_w(utmp,g)

    end do
  end do

  ! Setting boundary conditions
  bcs(:) = 2000

  ! Initializing solver
  call initialize(ns_solver,grid,dt,t_final,g,R,w0,winfty,bcs,"none",.true.,niter,tol,cfl)

  ! Solving the problem
  call solve_steady(ns_solver)

  ! Computing the L2 norm of the error

  ! Writing results to file

end program ns_validation
