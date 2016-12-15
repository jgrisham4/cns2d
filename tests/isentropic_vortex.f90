program isentropic_vortex
  use mesh_class,   only : mesh,read_from_file,preprocess
  use euler_solver, only : solver,initialize,solve_feuler,solve_rk4
  implicit none
  double precision :: a,rho_inf,p_inf,u_inf,v_inf,x0,y0,K,xb,yb,rb,temp,g
  double precision :: rgas,t_inf,a_inf,winfty(4)
  double precision, allocatable :: w0(:,:,:)
  double precision, parameter :: pi = 4.0d0*atan(1.0d0)
  integer :: aerr,i,j

  ! Creating mesh and solver objects
  type(mesh)   :: grid
  type(solver) :: esolver

  ! Reading mesh
  call read_from_file(grid,"iv800.cgns")
  !call read_from_file(grid,"iv200.cgns")

  ! Preprocessing mesh
  call preprocess(grid)

  ! Allocating memory for the initial field
  allocate(w0(grid%nelemi,grid%nelemj,4),stat=aerr)
  if (aerr.ne.0) then
    print *, "Can't allocate memory for w0 in isentropic vortex."
    stop
  end if

  ! Setting exact solution
  rgas      = 287.0d0
  g         = 1.4d0
  a         = 1.0d0
  rho_inf   = 1.0d0
  p_inf     = 1.0d0/1.4d0
  u_inf     = 2.0d0
  v_inf     = 2.0d0
  x0        = -10.0d0
  y0        = -10.0d0
  K         = 5.0d0
  t_inf     = p_inf/(rho_inf*rgas)
  a_inf     = sqrt(g*rgas*t_inf)
  winfty(1) = rho_inf
  winfty(2) = u_inf
  winfty(3) = v_inf
  winfty(4) = p_inf
  do j=1,grid%nelemj
    do i=1,grid%nelemi

      ! computing bar quantities
      xb = grid%elem(i,j)%xc - x0
      yb = grid%elem(i,j)%yc - y0
      rb = sqrt(xb**2 + yb**2)

      ! Assigning initial field
      temp = t_inf*(1.0d0 - K**2*(g-1.0d0)/(8.0d0*a*pi**2*a_inf**2)*exp(a*(1.0d0-rb**2)))
      w0(i,j,1) = rho_inf*(temp/t_inf)**(1.0d0/(g-1.0d0))
      w0(i,j,2) = u_inf - K/(2.0d0*pi)*yb*exp(a*(1.0d0-rb**2)/2.0d0)
      w0(i,j,3) = v_inf + K/(2.0d0*pi)*xb*exp(a*(1.0d0-rb**2)/2.0d0)
      w0(i,j,4) = p_inf*(temp/t_inf)**(g/(g-1.0d0))

    end do
  end do

  ! Initializing solver
  !call initialize(esolver,grid,0.001d0,10.0d0,1.4d0,w0,winfty)
  call initialize(esolver,grid,0.01d0,10.0d0,1.4d0,w0,winfty)

  ! Solving problem
  !call solve_feuler(esolver,1000)
  call solve_rk4(esolver,100)

end program isentropic_vortex
