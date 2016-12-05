program test_mesh
  use mesh_class, only : mesh, read_from_file, write_to_tec, preprocess
  implicit none
  character (len=30) :: mesh_file,tec_file

  ! Creating mesh object
  type(mesh) :: grid

  ! Defining inputs
  !mesh_file = "simple.cgns"
  !tec_file  = "simple.tec"
  mesh_file = "square.cgns"
  tec_file  = "square.tec"
  
  ! Reading mesh
  call read_from_file(grid,mesh_file)

  ! Preprocessing the mesh
  call preprocess(grid)

  ! Writing mesh to file
  call write_to_tec(grid,tec_file)

end program 
