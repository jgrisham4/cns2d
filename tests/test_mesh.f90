program test_mesh
  use class_mesh, only : mesh, read_from_file, write_to_tec
  implicit none
  character (len=30) :: mesh_file,tec_file

  ! Creating mesh object
  type(mesh) :: grid

  ! Defining inputs
  mesh_file = "simple.cgns"
  tec_file  = "simple.tec"
  
  ! Reading mesh
  call read_from_file(grid,mesh_file)

  ! Writing mesh to file
  call write_to_tec(grid,tec_file)

end program 
