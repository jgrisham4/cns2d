!===============================================================================
! This is the main program that reads inputs from the namelist cns2d.nml.
! It calls the necessary subroutines to run the solver.
!
! Author: James Grisham
! Date: 01/13/2017
!===============================================================================

program precns2d
  use mesh_class,   only : mesh,read_from_file,preprocess,compute_max_timesteps_inv,compute_max_timesteps_visc
  use solver_class, only : solver,initialize
  use temporal,     only : solve_feuler,solve_rk4,solve_steady
  implicit none
  logical                       :: viscous_terms
  character (len=30)            :: mesh_file,method,limiter
  double precision              :: u_inf,v_inf,p_inf,rho_inf
  double precision              :: g,C1,S,final_time,time_step,R,cfl,tol,dt_max_global,k
  integer                       :: bcids(4),write_freq,niter,niterfo
  type(mesh)                    :: grid
  type(solver)                  :: solv
  double precision              :: winfty(4)
  double precision, allocatable :: w0(:,:,:),dt(:,:)
  integer                       :: aer,i,j

  ! Namelists
  namelist /governing_equations/viscous_terms
  namelist /mesh_inputs/mesh_file
  namelist /freestream_properties/u_inf,v_inf,p_inf,rho_inf
  namelist /gas_properties/g,C1,S,R
  namelist /time_advancement/method,final_time,time_step,niter,tol,cfl,niterfo
  namelist /boundary_conditions/bcids
  namelist /slope_limiter/limiter,k
  namelist /output/write_freq

  ! Opening file
  open(unit=101,file="cns2d.nml",status="old")
  read(101,governing_equations)
  read(101,mesh_inputs)
  read(101,freestream_properties)
  read(101,gas_properties)
  read(101,time_advancement)
  read(101,boundary_conditions)
  read(101,slope_limiter)
  read(101,output)
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
  write(*,'(a,f12.5,a)') "u_inf  : ", u_inf, " m/s"
  write(*,'(a,f12.5,a)') "v_inf  : ", v_inf, " m/s"
  write(*,'(a,es12.5,a)') "p_inf  : ", p_inf, " Pa"
  write(*,'(a,f12.5,a/)') "rho_inf: ", rho_inf, " kg/m^3"
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Gas properties"
  write(*,'(a)') "-------------------------------"
  write(*,'(a,f5.2/)') "gamma : ", g
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Time advancement"
  write(*,'(a)') "-------------------------------"
  write(*,'(2a)') "method    : ", method
  write(*,'(a,f5.4,a)') "Final time: ", final_time, " sec"
  write(*,'(a,es12.5,a/)') "Time step : ", time_step, " sec"
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Boundary conditions"
  write(*,'(a)') "-------------------------------"
  write(*,'(a,i4)') "BC 1: ", bcids(1)
  write(*,'(a,i4)') "BC 2: ", bcids(2)
  write(*,'(a,i4)') "BC 3: ", bcids(3)
  write(*,'(a,i4/)') "BC 4: ", bcids(4)
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Slope limiter"
  write(*,'(a)') "-------------------------------"
  write(*,'(a,a/)') "limiter: ", limiter
  write(*,'(a)') "-------------------------------"
  write(*,'(a)') "Output options"
  write(*,'(a)') "-------------------------------"
  write(*,'(a,i9/)') "write_frequency: ", write_freq

  ! Setting freestream quantities
  winfty(1) = rho_inf
  winfty(2) = u_inf
  winfty(3) = v_inf
  winfty(4) = p_inf

  ! Reading mesh
  call read_from_file(grid,mesh_file)

  ! Preprocessing mesh
  call preprocess(grid)

  ! Need to add option to read restart from CGNS file
  ! Allocating memory for initial guess for solver
  allocate(w0(grid%nelemi,grid%nelemj,4),stat=aer)
  if (aer.ne.0) then
    print *, "Error: Can't allocate memory for w0 in precns2d."
    stop
  end if

  ! Filling in freestream values as initial guess
  do j=1,grid%nelemj
    do i=1,grid%nelemi
      w0(i,j,:) = winfty
    end do
  end do

  ! Initializing solver
  call initialize(solv,grid,time_step,final_time,g,R,w0,winfty,bcids,limiter,viscous_terms,niter,niterfo,tol,cfl,k)

  ! Computing the max timestep possible
  if (solv%is_visc.eq..true.) then
    call compute_max_timesteps_visc(solv%grid,solv%g,solv%R,cfl)
  else
    call compute_max_timesteps_inv(solv%grid,solv%g,solv%R,cfl)
  end if

  ! Finding the minimum timestep
  allocate(dt(grid%nelemi,grid%nelemj))
  do j=1,grid%nelemj
    do i=1,grid%nelemi
      dt(i,j) = solv%grid%elem(i,j)%dt_max
    end do
  end do
  dt_max_global = minval(dt)
  write(*,'(2(a,es15.8))') "Max dt for CFL = ", cfl, " is ", dt_max_global

end program precns2d
