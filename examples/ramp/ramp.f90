program isentropic_vortex
  use mesh_class,   only : mesh,read_from_file,preprocess
  use solvers, only : solver,initialize,solve_feuler,solve_rk4
  implicit none
  double precision :: a,rho_inf,p_inf,u_inf,v_inf,x0,y0,K,xb,yb,rb,temp,g
  double precision :: rgas,t_inf,a_inf,winfty(4),beta_r,y_shock
  double precision :: rho2,u2,v2,p2
  double precision :: dt,tfinal
  double precision, allocatable :: w0(:,:,:)
  double precision, parameter :: pi = 4.0d0*atan(1.0d0)
  integer :: aerr,i,j,bcs(4)


  ! Creating mesh and solver objects
  type(mesh)   :: grid
  type(solver) :: esolver

  ! Reading mesh
  call read_from_file(grid,"ramp200.cgns")

  ! Preprocessing mesh
  call preprocess(grid)

  ! Allocating memory for the initial field
  allocate(w0(grid%nelemi,grid%nelemj,4),stat=aerr)
  if (aerr.ne.0) then
    print *, "Can't allocate memory for w0."
    stop
  end if

  ! Time step and final time
  dt = 5.0e-10
  tfinal = 1.0e-5

  ! Setting initial guess
  !beta_r    = 29.314*pi/180.0d0
  beta_r    = 35.0d0*pi/180.0d0
  rgas      = 287.0d0
  g         = 1.4d0
  a         = 1.0d0
  rho_inf   = 0.1226d0
  p_inf     = 1.01325e4
  u_inf     = 670.0121d0
  v_inf     = -118.1412d0
  t_inf     = p_inf/(rho_inf*rgas)
  a_inf     = sqrt(g*rgas*t_inf)
  winfty(1) = rho_inf
  winfty(2) = u_inf
  winfty(3) = v_inf
  winfty(4) = p_inf
  rho2      = 0.1788d0
  u2        = 603.6766d0
  v2        = 0.0d0
  p2        = 1.729191e4
  do j=1,grid%nelemj
    do i=1,grid%nelemi

      ! Assigning initial field
      y_shock = grid%elem(i,j)%xc*tan(beta_r)

      if (grid%elem(i,j)%yc.ge.y_shock) then
        w0(i,j,1) = rho_inf
        w0(i,j,2) = u_inf
        w0(i,j,3) = v_inf
        w0(i,j,4) = p_inf
      else
        w0(i,j,1) = rho2
        w0(i,j,2) = u2
        w0(i,j,3) = v2
        w0(i,j,4) = p2
      end if

    end do
  end do

  ! Setting boundary conditions
  bcs(1) = 1003
  bcs(2) = 1001
  bcs(3) = 1000
  bcs(4) = 1000

  ! Initializing solver
  call initialize(esolver,grid,dt,tfinal,1.4d0,w0,winfty,bcs,"barth")
  !call initialize(esolver,grid,dt,tfinal,1.4d0,w0,winfty,bcs,"none")

  ! Solving problem
  call solve_feuler(esolver,1000)
  !call solve_rk4(esolver,500)

end program isentropic_vortex
