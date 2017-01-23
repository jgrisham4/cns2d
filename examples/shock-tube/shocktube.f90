program shocktube
  use utils,        only : u_to_w
  use mesh_class,   only : mesh,read_from_file,preprocess
  use solver_class, only : solver,initialize
  use temporal,     only : solve_feuler,solve_rk4
  implicit none
  double precision :: a,rho_inf,p_inf,u_inf,v_inf,x0,y0,K,xb,yb,rb,temp,g
  double precision :: rgas,t_inf,a_inf,winfty(4)
  double precision :: rho2,u2,v2,p2
  double precision :: dt,tfinal
  double precision, allocatable :: w0(:,:,:)
  double precision, parameter :: pi = 4.0d0*atan(1.0d0)
  integer :: aerr,i,j,bcs(4)

  ! Creating mesh and solver objects
  type(mesh)   :: grid
  type(solver) :: esolver

  ! Reading mesh
  call read_from_file(grid,"tube.cgns")

  ! Preprocessing mesh
  call preprocess(grid)

  ! Allocating memory for the initial field
  allocate(w0(grid%nelemi,grid%nelemj,4),stat=aerr)
  if (aerr.ne.0) then
    print *, "Can't allocate memory for w0."
    stop
  end if

  ! Time step and final time
  !dt = 1.0e-4
  dt = 5.0e-4
  tfinal = 0.2d0

  ! Setting initial guess
  winfty(1) = 0.0d0
  winfty(2) = 0.0d0
  winfty(3) = 0.0d0
  winfty(4) = 0.0d0
  do j=1,grid%nelemj
    do i=1,grid%nelemi

      if (grid%elem(i,j)%xc.le.0.5d0) then
        w0(i,j,1) = 1.0d0
        w0(i,j,2) = 0.0d0
        w0(i,j,3) = 0.0d0
        w0(i,j,4) = 1.0d0
      else
        w0(i,j,1) = 0.125d0
        w0(i,j,2) = 0.0d0
        w0(i,j,3) = 0.0d0
        w0(i,j,4) = 0.1d0
      end if

    end do
  end do

  ! Setting boundary conditions
  ! Weakly enforced slip wall
  bcs(1) = 1002
  bcs(2) = 1002
  bcs(3) = 1002
  bcs(4) = 1002

  ! Strongly enforced slip wall
  !bcs(1) = 1003
  !bcs(2) = 1003
  !bcs(3) = 1003
  !bcs(4) = 1003

  ! Initializing solver
  !call initialize(esolver,grid,dt,tfinal,1.4d0,287.0d0,w0,winfty,bcs,"barth",.false.,0,0,0.0d0,0.0d0,1.0d0)
  call initialize(esolver,grid,dt,tfinal,1.4d0,287.0d0,w0,winfty,bcs,"venkat",.false.,0,0,0.0d0,0.0d0,1.0d0)

  ! Solving problem
  !call solve_feuler(esolver,100)
  call solve_rk4(esolver,100)

  ! Writing out a slice to file
  j = 10
  !open(3,file="results2d.dat")
  open(3,file="./comparison/results2d-venkat.dat")
  do i=1,grid%nelemi
    esolver%grid%elem(i,j)%w = u_to_w(esolver%grid%elem(i,j)%u,esolver%g)
    write(3,'(es13.6,a,4(es13.6,1x))') esolver%grid%elem(i,j)%xc, " ", esolver%grid%elem(i,j)%w
  end do
  close(3)

end program shocktube
