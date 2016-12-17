program main
  implicit none
  logical            :: viscous_terms
  character (len=30) :: mesh_file,method
  double precision   :: u_inf,v_inf,p_inf,rho_inf
  double precision   :: g,C1,S,final_time,time_step

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
  write(*,*) "========================================"
  write(*,*) "Namelist information: "
  write(*,*) "========================================\n"
  write(*,*) "----- Governing equations -----"
  write(*,'(a)') "viscous terms: ", viscous_terms
  write(*,*) "----- Mesh inputs -----"
  write(*,'(a)') "mesh file    : ", mesh_file
  write(*,*) "----- Freestream properties -----"
  write(*,'(a,f5.2,a)') "u_inf        : ", u_inf, " m/s"
  write(*,'(a,f5.2,a)') "v_inf        : ", v_inf, " m/s"
  write(*,'(a,f5.2,a)') "p_inf        : ", p_inf, " Pa"
  write(*,'(a,f5.2,a)') "rho_inf      : ", rho_inf, " kg/m^3"
  write(*,*) "----- Gas properties -----"
  write(*,'(a,f5.2)') "gamma        : ", g
  write(*,*) "----- Time advancement -----"
  write(*,*) "method : ", method
  write(*,'(a,f5.2,a)') "Final time: ", final_time, " sec"
  write(*,'(a,f5.2,a)') "Time step : ", time_step, " sec"
  write(*,*) "----- Boundary conditions -----"
  write(*,'(a,i4)') "BC 1: ", bcids(1)
  write(*,'(a,i4)') "BC 2: ", bcids(2)
  write(*,'(a,i4)') "BC 3: ", bcids(3)
  write(*,'(a,i4)') "BC 4: ", bcids(4)


end program main
