program cns2d
  implicit none

  ! Modules
  use mesh_class, only: mesh,read_from_file,preprocess
  use solvers,    only: solver, initialize,solve_feuler,solve_rk4

  ! Namelist variables
  logical            :: viscous_terms
  character (len=30) :: mesh_file,method
  double precision   :: u_inf,v_inf,p_inf,rho_inf
  double precision   :: g,C1,S,final_time,time_step
  integer            :: bcids(4)

  ! Local variables
  type(mesh)   :: grid
  type(solver) :: s

  ! Namelists
  namelist /governing_equations/viscous_terms
  namelist /mesh/mesh_file
  namelist /freestream_properties/u_inf,v_inf,p_inf,rho_inf
  namelist /gas_properties/g,C1,S
  namelist /time_advancement/method,final_time,time_step
  namelist /boundary_conditions/bcids

  ! Opening file
  open(unit=101,file="cns2d.nml",status="old")
  read(101,governing_equations)
  read(101,mesh)
  read(101,freestream_properties)
  read(101,gas_properties)
  read(101,time_advancement)
  read(101,boundary_conditions)
  close(101)

  ! Printing some info
  write(*,'(/a)') "========================================"
  write(*,'(a)') "Namelist information: "
  write(*,'(a/)') "========================================"
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Governing equations"
  write(*,'(a)') "-------------------------------"
  write(*,'(a,l1/)') "viscous terms: ", viscous_terms
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Mesh inputs"
  write(*,'(a)') "-------------------------------"
  write(*,'(a,a/)') "mesh file: ", mesh_file
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Freestream properties"
  write(*,'(a)') "-------------------------------"
  write(*,'(a,f5.2,a)') "u_inf  : ", u_inf, " m/s"
  write(*,'(a,f5.2,a)') "v_inf  : ", v_inf, " m/s"
  write(*,'(a,f5.2,a)') "p_inf  : ", p_inf, " Pa"
  write(*,'(a,f5.2,a/)') "rho_inf: ", rho_inf, " kg/m^3"
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Gas properties"
  write(*,'(a)') "-------------------------------"
  write(*,'(a,f5.2/)') "gamma : ", g
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Time advancement"
  write(*,'(a)') "-------------------------------"
  write(*,'(2a)') "method    : ", method
  write(*,'(a,f5.2,a)') "Final time: ", final_time, " sec"
  write(*,'(a,f5.2,a/)') "Time step : ", time_step, " sec"
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Boundary conditions"
  write(*,'(a)') "-------------------------------"
  write(*,'(a,i4)') "BC 1: ", bcids(1)
  write(*,'(a,i4)') "BC 2: ", bcids(2)
  write(*,'(a,i4)') "BC 3: ", bcids(3)
  write(*,'(a,i4/)') "BC 4: ", bcids(4)

  ! Reading mesh
  call read_from_file(grid,mesh_file)

  ! Preprocessing mesh
  call preprocess(grid)

  ! Initializing solver
  call initialize(s,grid,time_step,final_time,g,w0,winfty,bcids)

  ! Solving the problem
  if (method.eq."forward_euler") then
    call solve_feuler(s,write_freq)
  else if (method.eq."rk4") then
    call solve_rk4(s,write_freq)
  else
    print *, "Error: Time advancement method ", method, " not recognized."
    stop
  end if

end program cns2d
