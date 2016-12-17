program main
  implicit none

  ! Namelists
  namelist /governing_equations/viscous_terms
  namelist /mesh/mesh_file
  namelist /freestream_properties/u_inf,v_inf,p_inf,rho_inf
  namelist /gas_properties/g,C1,S
  namelist /temporal_controls/time_advancement,final_time,time_step
  namelist /boundary_conditions/bcids


end program main
