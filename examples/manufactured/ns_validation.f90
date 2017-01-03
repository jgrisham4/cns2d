program ns_validation
  use mesh_class, only : mesh,read_from_file,preprocess
  use solvers,    only : solver,initialize,solve_feuler,solve_rk4
  implicit none
  double precision :: dt,t_final
  double precision, allocatable :: w0(:,:,:)
  integer :: aerr,i,j,bcs(4)

  ! Creating mesh and solver objects
  type(mesh)   :: grid
  type(solver) :: ns_solver

  ! Reading mesh
  call read_from_file(grid,"grid_41x31.cgns")

  ! Preprocessing grid
  call preprocess(grid)

  ! Allocating memory for the initial guess
  allocate(w0(grid%nelemi,grid%nelemj,4),stat=aerr)
  if (aerr.ne.0) then
    print *, "Can't allocate memory for w0."
    stop
  end if

  ! Setting time step and final time
  dt = 1.0e-5
  t_final = 0.2d0

end program ns_validation
